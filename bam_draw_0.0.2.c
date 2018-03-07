//gcc -g -O3 -Wall bam_draw_0.0.2.c -o bam_draw_0_0_2 ./samtools-0.1.19/libbam.a  -I./samtools-0.1.19 -L./samtools-0.1.19 -lpthread -lgd -lpng -lz -ljpeg -lfreetype -lm
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <err.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "sam.h"
#include "gd.h"
//#include "gdfontt.h"  //gdFontGetTiny()
//#include "gdfonts.h"//gdFontGetSmall()
//#include "gdfontmb.h" //gdFontGetMediumBold()
//#include "gdfontl.h"  //gdFontGetLarge()
#include "gdfontg.h"  //gdFontGetGiant()
#include "hash-master/hash.h"

#define MAX_DEPTH 500

#define LEFT 50
#define TOP 50
#define EACH_READ_HEIGHT 5
#define SCALE_WIDTH 16000

#define MAX_GENE_NUM 30000
#define MAX_CHR 64

typedef struct _data_buf_ {
	samfile_t *fp;
	bam_index_t *idx;
	unsigned int max_depth;
	struct list **my_list ;
	struct list **Stack ;
}DataBuf;

typedef struct _region_ {
	char *chromosome;
	int32_t pos;
	int32_t l_qseq;
	int32_t qlen;
	uint32_t end;
	uint16_t n_cigar;
	uint32_t *cigar;
	char *qname;
	char *qseq;
	char strand;
}Region;

typedef struct _gtfLine_
{
	char *seqname;
	char *feature;
	unsigned int start;
	unsigned int end ;
	char strand;
	char frame;
	char *attributes;
}GtfLine;

typedef struct _Gene_
{
	unsigned int start;
	unsigned int end ;
	struct list *gtflist;
}Gene;

typedef struct _Chr_
{
	unsigned int gene_count;
	Gene *chr_gene;
}Chr;

typedef struct _GeneCluster_
{
	unsigned int start;
	unsigned int end ;
	struct list *GeneList;
}GeneCluster;

typedef struct _ChrCluster_
{
	unsigned int Cluster_count;
	GeneCluster *chr_gene_cluster;
}ChrCluster;


struct globalArgs_t {
	char **infiles;
	unsigned short numInfiles;
	int overlap;
	int window;
	char *region;
	char *outfile;
	char *gtf;
	int strand;
} globalArgs;

void display_usage(char * argv[]);
static inline char *cal_GC(const uint8_t *seq,int32_t l_qseq,unsigned short *n_GC) ;
static inline int fetch_func(const bam1_t *b, void *data);
static inline void fill_item(const bam1_t *b,const bam1_core_t *c,uint32_t *cigar,char *qseq,Region *r,DataBuf *databuf);
static inline void parse_my_list(DataBuf *databuf,unsigned int j);
static inline unsigned int Stack_list(struct list *LIST,struct list **Stack);
gdImagePtr creat_png(DataBuf *databuf,unsigned short numInfiles,int *color,unsigned int *max_transcript_num);
static inline void draw_region(struct list **Stack,gdImagePtr im,unsigned int *current_top,unsigned int max_depth,int beg,int end,int color[12],const char *filename);
void write_png(gdImagePtr im,const char *chromosome,int beg,int end,const char *filename);
static inline void Region_delete_int(void *DATA);
static inline long long usec(void);

/*load GTF to a hash table*/
static inline void assign_gene_start_end(Gene *gene);
struct hash *load_gtf(const char *infile);
static inline void data_delete_listT(void *DATA) ;
static inline void data_delete_hashT(void *DATA) ;
static inline void print_GtfList(struct list *gtfline_list);
static inline void print_gtf(struct hash *hash_gtf);
static inline int cmp_gene(const void *a, const void *b);

/*cluster GTF*/
static inline int overlap(unsigned int S1,unsigned int E1,unsigned int S2,unsigned int E2);
struct hash *cluster_gtf(struct hash *hash_gtf);
static inline void data_delete_list_GeneT(void *DATA);
static inline void data_delete_hash_clusterT(void *DATA);

/*binary search GeneList*/
static inline void print_GeneList(struct list *GeneList);
static inline int test_bsearch(struct hash *Clusterd_GTF,DataBuf *databuf,int ref,int beg,int end);
static inline int bsearch_04(ChrCluster *chrcluster,int x);

static inline int bsearch_draw(ChrCluster *chr_culster,int beg,int end,unsigned int *max_transcript_num);
static inline void draw_GeneList_rectangle(unsigned short transcript_num,gdImagePtr im,unsigned int *current_top,int beg,int end);
static inline void draw_GeneList(GeneCluster *gene_cluster,unsigned short transcript_num,gdImagePtr im,unsigned int *current_top,int beg,int end,int color[12]);

static inline void draw_all_genes(struct hash *Clusterd_GTF,DataBuf *databuf,const char *current_dir,char *work_dir);
static inline void draw_a_region(struct hash *Clusterd_GTF,DataBuf *databuf);

long long usec(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(8092*sizeof(char));
	const char* usage=
"\nCopyright (c) 2012-2013\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Usage: %s [-w window_size] [-o overlap_size] [-g gtf_file][-r chr1:1-2000000] [-s 0] [-O OUTFILE] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for drawing a region(less than 500kb)'s reads distribution with samtools indexed bam files.\n" \
"Example1:\n  %s -w 20000 -g ~/data/mm10/mm10.gtf -r chr3:14500000-14570000 -O out /ifs5/ST_ANNO/USER/xuxiong/data/accepted_hits.bam\n" \
"Example2:\n  %s -w 20000 -g /home/xuxiong/data/hg19/hg19.gtf -r  chr1:150521898-150533413 -O out /share/work1/staff/xuxiong/twins/batch2/test_20131224/RT202FD_L2_I006.R1.clean.fastq.gz_1224/accepted_hits.bam\n" \
"\n" \
"   [-w window_size]   = window size. default is 20000.                 [option]\n" \
"   [-o overlap_size]  = overlap size. default is 0.                    [option]\n" \
"   [-g gtf_file]      = gtf_infile                                     [required]\n" \
"   [-r]               = region, default is whole genome.(chr1:1-20000) [option]\n" \
"   [-s]               = bool variant,strand or not, default is 0.      [option]\n" \
"   [-O OUTPUT_FILE]   = OUTPUT file. default is 'out'                  [option]\n" \
"   [-h]               = This helpful help screen.                      [option]\n" \
"   Infiles            = bam format input file(s),at least 1 bam file.  [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	exit(1);
}

char *cal_GC(const uint8_t *seq,int32_t l_qseq,unsigned short *n_GC) {
	char *temp_seq=(char *)malloc((l_qseq+1)*sizeof(char));
	memset(temp_seq,0,l_qseq+1);
	int32_t i=0;
	int temp=0;
	for (i=0;i<l_qseq ;i++ ) {
		temp=bam1_seqi(seq, i);
		switch (temp) {
			case 1:
				temp_seq[i]='A';
				break;
			case 2:
				temp_seq[i]='C';
				break;
			case 4:
				temp_seq[i]='G';
				break;
			case 8:
				temp_seq[i]='T';
				break;
			case 15:
				temp_seq[i]='N';
				break;
			default:
				break;
		}
		if (temp == 2 || temp == 4) (*n_GC)++;
	}
//	fprintf(stderr,"%s %d\n",temp_seq,*n_GC);
	return temp_seq;
}

int fetch_func(const bam1_t *b, void *data) {
	DataBuf *databuf=(DataBuf *)data;
	uint32_t *cigar = bam1_cigar(b);
	const bam1_core_t *c = &b->core;
	if (c->flag&BAM_FUNMAP) return 0;
//	fprintf(stdout,"%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%c\n",databuf->fp->header->target_name[c->tid],c->pos, bam_calend(c,cigar) ,c->pos + l,l,bam_cigar2qlen(c,cigar), bam1_qname(b), c->qual, c->bin,c->mtid,c->mpos,c->isize,(c->flag&BAM_FREVERSE)? '-' : '+');
	unsigned short current_GC =0;
	char *qseq=cal_GC(bam1_seq(b),c->l_qseq,&current_GC);
	Region *R=(Region *)malloc(sizeof(Region));
	fill_item(b,c,cigar,qseq,R,databuf);
//	fprintf(stdout,"%s\n",((Region *)list_tail(*(databuf->my_list+c->tid)))->qname);
	return 0;
}

void fill_item(const bam1_t *b,const bam1_core_t *c,uint32_t *cigar,char *qseq,Region *r,DataBuf *databuf){
	r->chromosome=strdup(databuf->fp->header->target_name[c->tid]);
	r->pos=c->pos;
	r->l_qseq=c->l_qseq;
	r->qlen=bam_cigar2qlen(c,cigar);
	r->end=bam_calend(c,cigar);
	r->cigar=malloc(sizeof(uint32_t)*c->n_cigar);
	memcpy(r->cigar,cigar,c->n_cigar*sizeof(uint32_t));
	r->n_cigar=c->n_cigar;
	r->qname=strdup(bam1_qname(b));
	r->qseq=qseq;
	r->strand=(c->flag&BAM_FREVERSE)? '-' : '+';
	list_push_back(*(databuf->my_list+c->tid), r);
}

void parse_my_list(DataBuf *databuf,unsigned int j){
	databuf->max_depth=0;
	if (!list_empty(*(databuf->my_list+j))) {
		databuf->Stack=(struct list **)calloc(MAX_DEPTH,sizeof(struct list *));
		databuf->max_depth = Stack_list(*(databuf->my_list+j),databuf->Stack);
	}
	else{
		list_delete(*(databuf->my_list+j));
	}
}

unsigned int Stack_list(struct list *LIST,struct list **Stack){
	struct item *current, *next;
	Region *r;
	unsigned int i=0;
	if (!LIST) return 0;
	unsigned int end=0;
	/******************differentiate junction******************
	while (1) {
		end=0;
		current = LIST->head;
		*(Stack+i)=list_create(Region_delete_int);
		while (current) {
			next = current->next;
			r = (Region *)current->data;
			if (r->n_cigar>1) {
				if (r->pos>end) {
					list_push_back(*(Stack+i), r);
					end=r->end;
					if (current->next && current->prev) {
						current->next->prev = current->prev;
						current->prev->next = next;
					}
					else if (current->prev) {
						LIST->tail->prev->next = NULL;
						LIST->tail=current->prev;
					}
					else if (current->next){
						LIST->head->next->prev = NULL;
						LIST->head = next;
					}
					else{
						LIST->tail->prev=NULL;
						LIST->tail->next=NULL;
						LIST->head->prev = NULL;
						LIST->head->next = NULL;
						LIST->head=current=LIST->tail;
					}
					if (current) {
						LIST->size--;
					}
				}
			}
			current = next ;
		}
		if (list_empty(*(Stack+i))) {
			list_delete(*(Stack+i));
			break;
		}
		if (++i > MAX_DEPTH) return i;
	}
	fprintf(stderr,"i: %d\tlist_size: %d\n",i,list_size(LIST));
	**********************************************************/
	while (!list_empty(LIST)) {
		end=0;
		current = LIST->head;
		*(Stack+i)=list_create(Region_delete_int);
		while (current) {
			next = current->next;
			r = (Region *)current->data;
			if (r->pos>end) {
				list_push_back(*(Stack+i), r);
				end=r->end;
				if (current->next && current->prev) {
					current->next->prev = current->prev;
					current->prev->next = next;
				}
				else if (current->prev) {
					LIST->tail->prev->next = NULL;
					LIST->tail=current->prev;
				}
				else if (current->next){
					LIST->head->next->prev = NULL;
					LIST->head = next;
				}
				else{
					LIST->tail->prev=NULL;
					LIST->tail->next=NULL;
					LIST->head->prev = NULL;
					LIST->head->next = NULL;
					LIST->head=current=LIST->tail;
				}
				if (current) {
					free(current);
					LIST->size--;
				}
			}
			current = next ;
		}
		if (++i > MAX_DEPTH) break;
	}
//	fprintf(stderr,"i: %d\tlist_size: %d\nDone stack list\n",i,list_size(LIST));
	if (!list_empty(LIST)) {
		list_delete(LIST);
	}else{
		free(LIST);
	}
	return i;
}

gdImagePtr creat_png(DataBuf *databuf,unsigned short numInfiles,int *color,unsigned int *max_transcript_num) {
	float scale=1;
	unsigned int Height = TOP,i=0,each_transcript_height=10;
	for (i=0;i<numInfiles ;i++ ) {
		Height+=1.2*EACH_READ_HEIGHT*(databuf+i)->max_depth+((databuf+i)->max_depth>0?TOP:0);
//		fprintf(stderr,"max_depth: %d\n",(databuf+i)->max_depth);
	}
	Height+=(*max_transcript_num)*each_transcript_height+(*max_transcript_num>0?TOP:0);
	unsigned int Width = SCALE_WIDTH/scale + 2*LEFT;
	gdImagePtr im = gdImageCreate(Width,Height);
	int rgb[12][3]={
			{255,255,255},	/*****background_color******/
//			{70, 130, 180},
			{255, 140, 0},{160, 82, 45},{135, 206, 235},{107, 142, 35},{106, 90, 205},
			{119, 136, 153},{218, 165, 32},{178, 34, 34},{255, 0, 255},{0, 255, 255},{0, 255, 0},
		};
	for (i=0;i<12 ;i++ ) {
		color[i] = gdImageColorExact(im,rgb[i][0], rgb[i][1],rgb[i][2]);
		if (color[i] != (-1)) gdImageColorDeallocate(im, color[i]);
		color[i] = gdImageColorAllocate(im,rgb[i][0], rgb[i][1],rgb[i][2]);
	}
	return im;
}

void draw_GeneList_rectangle(unsigned short transcript_num,gdImagePtr im,unsigned int *current_top,int beg,int end){
	float scale=1;
	unsigned int each_transcript_height=10;
	unsigned int Height = transcript_num*(each_transcript_height+5);
	unsigned int Width = SCALE_WIDTH/scale + 2*LEFT;

	int im_black = gdImageColorAllocate(im,0,0,0);
	gdImageSetThickness(im, 1);
	gdImageRectangle(im, LEFT,*current_top,Width-LEFT,*current_top+Height, im_black);
	(*current_top)+=Height+TOP;
}

void draw_GeneList(GeneCluster *gene_cluster,unsigned short transcript_num,gdImagePtr im,unsigned int *current_top,int beg,int end,int color[12]) {
	struct item *current, *next;

	float scale=1;
	unsigned int each_transcript_height=10;
	float ratio= (float)SCALE_WIDTH/(end-beg)/scale;

	float x=LEFT;
	unsigned int y=*current_top;
//	gdImageString(im,gdFontGetGiant(),x,y-TOP/2,(unsigned char *)filename,im_black);
	current = (gene_cluster->GeneList)->head;
	while (current) {
		y+=5;
		x=LEFT;
		next=current->next;
		Gene *current_gene=(Gene *)current->data;
//		fprintf(stderr,"%u\t%u\n",current_gene->start,current_gene->end);
		int current_color=(((GtfLine *)current_gene->gtflist->head->data)->strand == '+' ? color[1] : color[10]);
		gdImageSetThickness(im, 1);
		gdImageLine(im, (int)(x+(current_gene->start-beg)*ratio),y, (int)(x+(current_gene->end-beg)*ratio),y,current_color) ;
		struct item *current2,*next2;
		current2=current_gene->gtflist->head;
		while (current2) {
			next2=current2->next;
			GtfLine *gtfline=(GtfLine *)current2->data;
			if (!strcmp(gtfline->feature,"exon")) {
				gdImageSetThickness(im, 10);
				gdImageLine(im,(int)(x+(gtfline->start-beg)*ratio),y,(int)(x+(gtfline->end-beg)*ratio),y,current_color) ;
//				fprintf(stderr,"%s\t%s\t%u\t%u\t%c\t%c\t%s\n",gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
//				gdImageRectangle(im, (int)(x+(gtfline->start-beg)*ratio),y,(int)(x+(gtfline->end-beg)*ratio),y+each_transcript_height,current_color);
//				gdImageFilledRectangle(im,(int)(x+(gtfline->start-beg)*ratio),y,(int)(x+(gtfline->end-beg)*ratio),y+each_transcript_height,current_color);
			}
			current2=next2;
		}
		y+=each_transcript_height;
		current=next;
	}
}

void draw_region(struct list **Stack,gdImagePtr im,unsigned int *current_top,unsigned int max_depth,int beg,int end,int color[12],const char *filename) {
	unsigned int i=0,j=0,l=0;
	struct item *current, *next;
	Region *r;

	float scale=1;
	unsigned int Height = 1.2*max_depth*EACH_READ_HEIGHT;
	unsigned int Width = SCALE_WIDTH/scale + 2*LEFT;
	float ratio= (float)SCALE_WIDTH/(end-beg)/scale;
	
	int im_black = gdImageColorAllocate(im,0,0,0);
	gdImageSetThickness(im, 1);
	gdImageRectangle(im, LEFT,*current_top,Width-LEFT,*current_top+Height, im_black);

	float x=LEFT;
	unsigned int y=*current_top;
	gdImageString(im,gdFontGetGiant(),x,y-TOP/2,(unsigned char *)filename,im_black);
	
	for (i=0;i<max_depth ;i++ ) {
		x=LEFT;
		y+=EACH_READ_HEIGHT;
		if (list_empty(*(Stack+i))) break;
		current = (*(Stack+i))->head;
		while (current) {
			next = current->next;
			r = (Region *)current->data;
			if (current) {
				unsigned int temp_start= r->pos;
				for (j = 0,l = 0; j < r->n_cigar; ++j) {
					int styleDotted[2]={r->strand == '+' ? color[1] : color[10],gdTransparent};
					int op = r->cigar[j]&0xf;
					if ( op == BAM_CREF_SKIP){
						l = r->cigar[j]>>4;
						gdImageSetThickness(im, 1);
						gdImageSetStyle(im, styleDotted, 2);
					}
					else if (op == BAM_CDEL ){
						int styledel[2]={r->strand == '+' ? color[2] : color[9],r->strand == '+' ? color[2] : color[9]};
						l = r->cigar[j]>>4;
						gdImageSetThickness(im, 4);
						gdImageSetStyle(im, styledel, 1);
					}
					else if (op == BAM_CINS ){
						int styleins[2]={r->strand == '+' ? color[3] : color[8],r->strand == '+' ? color[3] : color[8]};
						l = r->cigar[j]>>4;
						gdImageSetThickness(im, 4);
						gdImageSetStyle(im, styleins, 1);
					}
					else if (op == BAM_CMATCH) {
						l = r->cigar[j]>>4;
						gdImageSetThickness(im, 4);
						gdImageSetStyle(im, styleDotted, 1);
					}
					gdImageLine(im, (int)(x+(temp_start-beg+1)*ratio), (int)y, (int)(x+(temp_start+l-beg+1)*ratio), (int)y,gdStyled) ;
					temp_start+=l;
				}
//				fprintf(stderr,"%d\t%d\t%d\t%d\t#%d#\n",i,r->pos,r->end,(*(Stack+i))->size,r->n_cigar);
				item_delete(*(Stack+i), current);
				(*(Stack+i))->size--;
			}
			current = next ;
		}
		free(*(Stack+i));
	}
	(*current_top)+=Height+TOP;
	free(Stack);
}

void write_png(gdImagePtr im,const char *chromosome,int beg,int end,const char *filename){
	char *output=(char *)malloc(1024*sizeof(char));
	memset(output,0,1024*sizeof(char));
	sprintf(output,"%s_%s_%d_%d.png",filename,chromosome,beg,end);
	FILE *reg_png=fopen(output,"wb");
	if(reg_png == NULL) err(1, "Failed to open file (%s)", output);
	gdImagePng(im,reg_png);
	gdImageDestroy(im);
	fclose(reg_png);
	free(output);
}

void Region_delete_int(void *DATA) {
	Region *R=(Region *)DATA;
	free(R->cigar);
	free(R->chromosome);
	free(R->qname);
	free(R->qseq);
	free(R);
}

void data_delete_hashT(void *DATA) {
	Gene *chr_gene=(Gene *)DATA;
	free(chr_gene);
}

void data_delete_listT(void *DATA) {
	GtfLine *gtfline=(GtfLine *)DATA;
	free(gtfline->seqname);
	free(gtfline->feature);
	free(gtfline->attributes);
	free(gtfline);
}

struct hash *load_gtf(const char *infile){
	FILE *f=fopen(infile,"r");
	if(f == NULL) {
		printf("Failed to open file (%s)", infile);
		exit(1);
	}
	char *buf=(char *)calloc(1024,sizeof(char));
	char *temp_seqname = "-";
	char *temp_attributes = "-";
	struct hash *hash_gtf = hash_create(data_delete_hashT, MAX_CHR);
	unsigned int i=0;
	Gene *temp_chr_gene=(Gene *)calloc(MAX_GENE_NUM,sizeof(Gene));
	while (1) {
		buf=fgets(buf,1024*sizeof(char),f);
		if (feof(f)) break;
		GtfLine *gtfline=(GtfLine *)malloc(sizeof(GtfLine));
		gtfline->seqname=(char *)calloc(1024,sizeof(char));
		gtfline->feature=(char *)calloc(1024,sizeof(char));
		gtfline->attributes=(char *)calloc(1024,sizeof(char));
		sscanf(buf, "%s\t%*[^\t]\t%s\t%u\t%u\t%*[^\t]\t%c\t%c\tgene_id \"%s",gtfline->seqname,gtfline->feature, &gtfline->start, &gtfline->end,&gtfline->strand,&gtfline->frame,gtfline->attributes);
		gtfline->seqname=realloc(gtfline->seqname,(strlen(gtfline->seqname)+1)*sizeof(char));
		gtfline->feature=realloc(gtfline->feature,(strlen(gtfline->feature)+1)*sizeof(char));
		gtfline->attributes=realloc(gtfline->attributes,(strlen(gtfline->attributes)+1)*sizeof(char));
//		fprintf(stderr,"%s\t%s\t%u\t%u\t%c\t%c\t%s\n",gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
	
		if (strcmp(temp_attributes,gtfline->attributes)) {
			if (strcmp(temp_attributes,"-")) {
				assign_gene_start_end(&temp_chr_gene[i]);
				i++;
			}
			if (strcmp(temp_seqname,gtfline->seqname)) {
				if (strcmp(temp_seqname,"-")) {
					Chr *CHR=(Chr *)malloc(sizeof(Chr));
					CHR->gene_count=i;
					CHR->chr_gene=(Gene *)calloc(CHR->gene_count,sizeof(Gene));
					memmove(CHR->chr_gene,temp_chr_gene,(CHR->gene_count)*sizeof(Gene));
					hash_insert(hash_gtf,temp_seqname,CHR);
//					fprintf(stderr,"%s: %d\n",temp_seqname,i);
				}
				i=0;
			}
			temp_chr_gene[i].gtflist=list_create(data_delete_listT);
			list_push_back(temp_chr_gene[i].gtflist,gtfline);
		}
		else{
			list_push_back(temp_chr_gene[i].gtflist,gtfline);
		}
		temp_seqname=gtfline->seqname;
		temp_attributes=gtfline->attributes;
	}

	Chr *CHR=(Chr *)malloc(sizeof(Chr));
	CHR->gene_count=++i;
	CHR->chr_gene=(Gene *)calloc(CHR->gene_count,sizeof(Gene));
	memmove(CHR->chr_gene,temp_chr_gene,(CHR->gene_count)*sizeof(Gene));
	hash_insert(hash_gtf,temp_seqname,CHR);
//	fprintf(stderr,"%s: %d",temp_seqname,i);
	fclose(f);
	free(temp_chr_gene);
	return hash_gtf;
}

void assign_gene_start_end(Gene *gene){
	struct item *current, *next,*prev;
	current = gene->gtflist->head;
	while (current) {
		next = current->next;
		GtfLine *temp_gtfline=(GtfLine *)current->data;
		if (!strcmp(temp_gtfline->feature,"exon")) {
			gene->start = temp_gtfline->start;
			break;
		}
		current = next;
	}
	current = gene->gtflist->tail;
	while (current) {
		prev = current->prev;
		GtfLine *temp_gtfline=(GtfLine *)current->data;
		if (!strcmp(temp_gtfline->feature,"exon")) {
			gene->end = temp_gtfline->end;
			break;
		}
		current = prev;
	}
}

void print_GtfList(struct list *gtfline_list){
	struct item *current, *next;
	current = gtfline_list->head;
	while (current) {
		next = current->next;
		GtfLine *gtfline=(GtfLine *)current->data;
		fprintf(stderr,"%s\t%s\t%u\t%u\t%c\t%c\t%s\n",gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
		current = next;
	}
}

void print_gtf(struct hash *hash_gtf){
	unsigned i=0,j=0;
	struct item *current, *next;
	fprintf(stdout,"%d\n",hash_gtf->key_count);
	for (i=0;i<hash_gtf->size ;i++ ) {
		struct list *list = (hash_gtf->table)[i];
//		fprintf(stdout,"bucket: %d\tlist_size: %d\n",i,list_size(list));
		current = list->head;
		while (current) {
			next = current->next;
			Chr *chr=(Chr *)((struct pair *)current->data)->value;
			fprintf(stdout,"%s: %u\n",(char *)((struct pair *)current->data)->key,chr->gene_count);
			qsort(chr->chr_gene,chr->gene_count,sizeof(Gene),cmp_gene);
			for (j=0;j< chr->gene_count;j++ ) {
				Gene current_gene=chr->chr_gene[j];
				fprintf(stderr,"%u\t%u\n",current_gene.start,current_gene.end);
				print_GtfList(current_gene.gtflist);
			}
			current = next;
		}
	}
}

int cmp_gene(const void *a, const void *b){
	Gene *c=(Gene *)a;
	Gene *d=(Gene *)b;
	if (c->start!=d->start) {
		return c->start-d->start;
	}
	else{
		return c->end-d->end;
	}
}

void data_delete_hash_clusterT(void *DATA) {
	ChrCluster *chrcluster=(ChrCluster *)DATA;
	free(chrcluster);
}

void data_delete_list_GeneT(void *DATA) {
	Gene *gene=(Gene *)DATA;
	list_delete(gene->gtflist);
	free(gene);
}

int overlap(unsigned int S1,unsigned int E1,unsigned int S2,unsigned int E2){
	return E1>=S2?1:0;
}

struct hash *cluster_gtf(struct hash *hash_gtf) {
	struct hash *hash_cluster_gtf=hash_create(data_delete_hash_clusterT, MAX_CHR);
	unsigned i=0,j=0;
	struct item *current, *next;
	for (i=0;i<hash_gtf->size ;i++ ) {
		struct list *list = (hash_gtf->table)[i];
		current = list->head;
		while (current) {
			next = current->next;
			Chr *chr=(Chr *)((struct pair *)current->data)->value;
//			fprintf(stdout,"%s: %u\n",(char *)((struct pair *)current->data)->key,chr->gene_count);
			qsort(chr->chr_gene,chr->gene_count,sizeof(Gene),cmp_gene);

			unsigned int cluster_index=0;
			GeneCluster *chr_gene_cluster=(GeneCluster *)calloc(MAX_GENE_NUM,sizeof(GeneCluster));
			chr_gene_cluster[cluster_index].start=0;
			chr_gene_cluster[cluster_index].end=0;

			for (j=0;j< chr->gene_count;j++ ) {
				if (!overlap(chr_gene_cluster[cluster_index].start,chr_gene_cluster[cluster_index].end,chr->chr_gene[j].start,chr->chr_gene[j].end)) {
					if (!list_empty(chr_gene_cluster[cluster_index].GeneList)) cluster_index++;
					chr_gene_cluster[cluster_index].GeneList=list_create(data_delete_list_GeneT);
					list_push_back(chr_gene_cluster[cluster_index].GeneList,&chr->chr_gene[j]);
					chr_gene_cluster[cluster_index].start=((Gene *)chr_gene_cluster[cluster_index].GeneList->head->data)->start;
				}
				else{
					list_push_back(chr_gene_cluster[cluster_index].GeneList,&chr->chr_gene[j]);
				}
				chr_gene_cluster[cluster_index].end=chr_gene_cluster[cluster_index].end>((Gene *)chr_gene_cluster[cluster_index].GeneList->tail->data)->end?chr_gene_cluster[cluster_index].end:((Gene *)chr_gene_cluster[cluster_index].GeneList->tail->data)->end;
			}
			
			ChrCluster *chr_culster=(ChrCluster *)malloc(sizeof(ChrCluster));
			chr_culster->Cluster_count=cluster_index+1;
//			fprintf(stderr,"%s\t%u\n",(char *)((struct pair *)current->data)->key,chr_culster->Cluster_count);
			chr_culster->chr_gene_cluster=(GeneCluster *)calloc(chr_culster->Cluster_count,sizeof(GeneCluster));
			memmove(chr_culster->chr_gene_cluster,chr_gene_cluster,chr_culster->Cluster_count*sizeof(GeneCluster));
			hash_insert(hash_cluster_gtf,((struct pair *)current->data)->key,chr_culster);
			free(chr_gene_cluster);
			current = next;
		}
	}
	return hash_cluster_gtf;
}

int bsearch_04(ChrCluster *chrcluster,int x) {
	int lo = 0;
	int hi = chrcluster->Cluster_count - 1;

	while(lo <= hi){
		int mid = (hi + lo) / 2;
//		int mid = (hi + lo) >> 1;
		if(x < chrcluster->chr_gene_cluster[mid].start)
			hi = mid - 1;
		else if(x > chrcluster->chr_gene_cluster[mid+1].start)
			lo = mid + 1;
		else if (x >=chrcluster->chr_gene_cluster[mid].start && x <=chrcluster->chr_gene_cluster[mid+1].start) {
			return mid;
		}
	}
	return -1;
}

void print_GeneList(struct list *GeneList){
	struct item *current, *next;
	if (!GeneList) return;
	fprintf(stderr,"GeneListSize: %d\n",list_size(GeneList));
	// Iterate through list and delete all items
	current = GeneList->head;
	while (current) {
		next = current->next;
		Gene *current_gene=(Gene *)current->data;
		fprintf(stderr,"Gene_start: %u\tGene_end: %u\n",current_gene->start,current_gene->end);
		print_GtfList(current_gene->gtflist);
		current = next;
	}
}

int test_bsearch(struct hash *Clusterd_GTF,DataBuf *databuf,int ref,int beg,int end) {
	ChrCluster *chr_culster=(ChrCluster *)hash_find(Clusterd_GTF, databuf->fp->header->target_name[ref]);
	int gene_cluster_index=bsearch_04(chr_culster,beg);
	if (gene_cluster_index>0) {
		fprintf(stderr,"gene_cluster_index: %d\tgene_cluster_start: %u\tgene_cluster_end: %u\n",gene_cluster_index,chr_culster->chr_gene_cluster[gene_cluster_index].start,chr_culster->chr_gene_cluster[gene_cluster_index].end);
		print_GeneList(chr_culster->chr_gene_cluster[gene_cluster_index].GeneList);
		while (chr_culster->chr_gene_cluster[gene_cluster_index].end <end) {
			print_GeneList(chr_culster->chr_gene_cluster[++gene_cluster_index].GeneList);
		}
	}
	return gene_cluster_index;
}

int bsearch_draw(ChrCluster *chr_culster,int beg,int end,unsigned int *max_transcript_num) {
	int gene_cluster_index=bsearch_04(chr_culster,beg);
	if (gene_cluster_index>0) {
		int i=gene_cluster_index;
		*max_transcript_num=list_size(chr_culster->chr_gene_cluster[i].GeneList);
		while (chr_culster->chr_gene_cluster[i].end <end) {
			int current_list_size = list_size(chr_culster->chr_gene_cluster[++i].GeneList);
			if (current_list_size>*max_transcript_num) {
				*max_transcript_num=current_list_size;
			}
		}
		fprintf(stderr,"max_transcript_num: %u\n",*max_transcript_num);
	}
	return gene_cluster_index;
}

void draw_a_region(struct hash *Clusterd_GTF,DataBuf *databuf){
	int ref, beg, end,i,boolen=0;
	int color[12];
	unsigned int current_top=50;
//	fprintf(stderr,"%s\t%d\t%d\n",databuf->fp->header->target_name[ref],beg,end);
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		ref=0;beg=0;end=0;
		bam_parse_region((databuf+i)->fp->header, globalArgs.region, &ref, &beg, &end);
		if (ref < 0) err(1, "bam2bed: Invalid region %s\n", globalArgs.region);
		*((databuf+i)->my_list+ref)=list_create(Region_delete_int);
		bam_fetch((databuf+i)->fp->x.bam, (databuf+i)->idx, ref, beg, end, (databuf+i), fetch_func);
		fprintf(stderr,"%s read count in this region: %d\n",*(globalArgs.infiles+i),list_size(*((databuf+i)->my_list+ref)));
		if (!list_empty(*((databuf+i)->my_list+ref))) boolen=1;
		parse_my_list(databuf+i,ref);
	}
	if (boolen) {
		unsigned int max_transcript_num=0;
//		int gene_cluster_index = test_bsearch(Clusterd_GTF,databuf,ref,beg,end);
//		fprintf(stderr,"Done binary search GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
		ChrCluster *chr_culster=(ChrCluster *)hash_find(Clusterd_GTF, databuf->fp->header->target_name[ref]);
		int gene_cluster_index=bsearch_draw(chr_culster,beg,end,&max_transcript_num);
		gdImagePtr im=creat_png(databuf,globalArgs.numInfiles,color,&max_transcript_num);
		if (gene_cluster_index>0) {
			draw_GeneList(&chr_culster->chr_gene_cluster[gene_cluster_index],max_transcript_num,im,&current_top,beg,end,color);
			while (chr_culster->chr_gene_cluster[gene_cluster_index].end <end) {
				draw_GeneList(&(chr_culster->chr_gene_cluster[++gene_cluster_index]),max_transcript_num,im,&current_top,beg,end,color);
			}
			draw_GeneList_rectangle(max_transcript_num,im,&current_top,beg,end);
		}
		for (i=0;i<globalArgs.numInfiles ;i++ ) {
			if ((databuf+i)->max_depth>0) {
				draw_region((databuf+i)->Stack,im,&current_top,(databuf+i)->max_depth,beg,end,color,*(globalArgs.infiles+i));
			}
		}
		write_png(im,(databuf+i-1)->fp->header->target_name[ref],++beg,end,globalArgs.outfile);
	}
}

void draw_all_genes(struct hash *Clusterd_GTF,DataBuf *databuf,const char *current_dir,char *work_dir){
	unsigned int i,j,m=0;
	int color[12];
	for (j=0;j<databuf->fp->header->n_targets ;j++ ) {
		sprintf(work_dir,"%s/%s",current_dir,databuf->fp->header->target_name[j]);
		mkdir(work_dir,0750);
		int status = chdir(work_dir);
		if (status) err(1,"entering directory error!\n");
		ChrCluster *chr_culster=(ChrCluster *)hash_find(Clusterd_GTF, databuf->fp->header->target_name[j]);
		
		if (chr_culster==NULL) continue;
		fprintf(stderr,"%s\t%u\n",databuf->fp->header->target_name[j],chr_culster->Cluster_count);
		for (m=0;m<chr_culster->Cluster_count ;m++ ) {
			int boolen=0;
			for (i=0;i<globalArgs.numInfiles ;i++ ) {
				*((databuf+i)->my_list+j)=list_create(Region_delete_int);
				bam_fetch((databuf+i)->fp->x.bam, (databuf+i)->idx,j,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,(databuf+i), fetch_func);
				fprintf(stderr,"%u\t%s read count in this region: %d\n",m+1,*(globalArgs.infiles+i),list_size(*((databuf+i)->my_list+j)));
				if (!list_empty(*((databuf+i)->my_list+j))) boolen=1;
				parse_my_list(databuf+i,j);
			}
			if (boolen) {
				unsigned int max_transcript_num=list_size(chr_culster->chr_gene_cluster[m].GeneList);
				unsigned int current_top=50;

				ChrCluster *chr_culster=(ChrCluster *)hash_find(Clusterd_GTF, databuf->fp->header->target_name[j]);
				gdImagePtr im=creat_png(databuf,globalArgs.numInfiles,color,&max_transcript_num);
				draw_GeneList(&chr_culster->chr_gene_cluster[m],max_transcript_num,im,&current_top,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,color);
				draw_GeneList_rectangle(max_transcript_num,im,&current_top,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end);
				
				for (i=0;i<globalArgs.numInfiles ;i++ ) {
					if ((databuf+i)->max_depth>0) {
						draw_region((databuf+i)->Stack,im,&current_top,(databuf+i)->max_depth,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,color,*(globalArgs.infiles+i));
					}
				}
				write_png(im,(databuf+i-1)->fp->header->target_name[j],chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,globalArgs.outfile);
			}
		}
	}
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.overlap=0;
	globalArgs.window=20000;
	globalArgs.region="-";
	globalArgs.outfile="out";
	globalArgs.strand=0;
	globalArgs.gtf="-";

	const char *optString = "o:w:g:r:s:O:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.overlap = atoi(optarg);
				break;
			case 'w':
				globalArgs.window = atoi(optarg);
				break;
			case 'r':
				globalArgs.region = optarg;
				break;
			case 's':
				globalArgs.strand = atoi(optarg);
				break;
			case 'g':
				globalArgs.gtf = optarg;
				break;
			case 'O':
				globalArgs.outfile = optarg;
			case '?':	/* fall-through is intentional */
				break;
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	globalArgs.infiles = argv + optind;
	globalArgs.numInfiles = argc - optind;

	long long begin;
	begin=usec();

	struct hash *GTF=load_gtf(globalArgs.gtf);
	fprintf(stderr,"Done load GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
//	print_gtf(GTF);
	struct hash *Clusterd_GTF=cluster_gtf(GTF);
	fprintf(stderr,"Done clusted GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));

	unsigned int i=0;
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			err(1, "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
		}
		(databuf+i)->my_list = (struct list **)calloc((databuf+i)->fp->header->n_targets,sizeof(struct list *));
		if (((databuf+i)->idx = bam_index_load(*(globalArgs.infiles+i))) == 0) {
			err(1, "bam2bed: BAM indexing file is not available.\n");
		}
//		fprintf(stderr,"$%d\t%d\n",j,(databuf+i)->fp->header->target_len[j]);
	}

	char *current_dir=(char *)calloc(128,sizeof(char));
	current_dir = getcwd(current_dir, 128*sizeof(char));
	char *work_dir=(char *)calloc(256,sizeof(char));
	
	if (strncmp(globalArgs.region,"-",1)==0) { /* if a region is not specified */
		draw_all_genes(Clusterd_GTF,databuf,current_dir,work_dir);
	}
	else{
		draw_a_region(Clusterd_GTF,databuf);
	}
	int status=chdir(current_dir);
	if (status) exit(0);
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		bam_index_destroy((databuf+i)->idx);
		free((databuf+i)->my_list);
		samclose((databuf+i)->fp);
	}
	free(databuf);
	free(current_dir);
	free(work_dir);

	fprintf(stderr,"Finished at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
