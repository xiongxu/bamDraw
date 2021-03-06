// gcc -g -O3 -Wall bam_draw_0.0.5.c -o bam_draw_0_0_5 -I.. -I./samtools-0.1.19 -L.. -L./samtools-0.1.19 $(pkg-config --cflags --libs cairo) -lpthread -lbam -lhashtbl
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <cairo.h>
#include <cairo-svg.h>
#include <cairo-pdf.h>
#include "sam.h"
#include "../hashtbl.h"
#include "../gtf.h"

#define MAX_DEPTH 200

#define LEFT 50
#define TOP 50
#define EACH_READ_HEIGHT 5
#define EACH_WIG_HEIGHT 5
#define SCALE_WIDTH 16000

typedef struct _data_buf_ {
	samfile_t *fp;
	bam_index_t *idx;
	unsigned int max_depth;
	double max_wig;
	struct list **my_list ;
	struct list **wig_list ;
	struct list **Stack ;
} DataBuf;

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
	char *MD;
} Region;

typedef struct _wig_
{
	int start;
	int end;
	double depth;
} Wig;

struct globalArgs_t {
	char **infiles;
	unsigned short numInfiles;
	int overlap;
	int depth;
	int pileup;
	char *region;
	char *outfile;
	char *gtf;
	int log;
} globalArgs;

static int rgb[12][3]={
			{255,255,255},
//			{70, 130, 180},
			{255, 140, 0},{160, 82, 45},{135, 206, 235},{107, 142, 35},{106, 90, 205},
			{119, 136, 153},{218, 165, 32},{178, 34, 34},{255, 0, 255},{0, 255, 255},{0, 255, 0},
		};

void display_usage(char * argv[]);
static inline long long usec(void);
static inline char *cal_GC(const uint8_t *seq,int32_t l_qseq,unsigned short *n_GC) ;
static inline int fetch_func(const bam1_t *b, void *data);
static inline void fill_item(const bam1_t *b,const bam1_core_t *c,uint32_t *cigar,char *qseq,Region *r,DataBuf *databuf);
static inline void parse_my_list(DataBuf *databuf,unsigned int j);
static inline unsigned int Stack_list(struct list *LIST,struct list **Stack);

static inline char *create_pdf(const char *chromosome,int beg,int end,const char *filename);
static inline void creat_pdf(DataBuf *databuf,unsigned short numInfiles,unsigned int *max_transcript_num,int beg,int end,int *Width,int *Height);
static inline void draw_GeneList_rectangle(unsigned short transcript_num,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,int beg,int end);
static inline void draw_GeneList(GeneCluster *gene_cluster,unsigned short transcript_num,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,int beg,int end);
static inline void draw_region(struct list **Stack,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,unsigned int max_depth,int beg,int end,const char *filename) ;
static inline void free_cario(cairo_t *cr,cairo_surface_t *surface,char *filename);

static inline void Region_delete_int(void *DATA);

double tity_sam2wig(struct list *LIST,struct list *wig_list);
void Wig_delete_int(void *DATA);
void draw_wig(struct list *wig_list,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,double max_wig,int beg,int end,const char *filename);

/*binary search GeneList*/
static inline void print_GeneList(struct list *GeneList);
static inline int test_bsearch(HASHTBL *Clusterd_GTF,DataBuf *databuf,int ref,int beg,int end);
static inline int bsearch_04(ChrCluster *chrcluster,int x);
static inline int bsearch_draw(ChrCluster *chr_culster,int beg,int end,unsigned int *max_transcript_num);

static inline void draw_all_genes(HASHTBL *Clusterd_GTF,DataBuf *databuf,const char *current_dir,char *work_dir);
static inline void draw_a_region(HASHTBL *Clusterd_GTF,DataBuf *databuf);

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
"Usage: %s [-D] [-P] [-o overlap_size] [-g gtf_file][-r chr1:1-2000000] [-s 0] [-O OUTFILE] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for drawing a region(less than 500kb)'s reads distribution with samtools indexed bam files.\n" \
"Example:\n  %s -D -P -g /home/xuxiong/data/hg19/hg19.gtf  -r chr1:762971-794826 -O out /share/work1/staff/xuxiong/twins/batch2/test_20131224/*_1224/accepted_hits.bam\n" \
"\n" \
"   [-D depth]         = draw bedgraph format depth.                     [option]\n" \
"   [-P reads pileup]  = draw reads pileup.                              [option]\n" \
"   [-o overlap_size]  = overlap size. default is 0.                     [option]\n" \
"   [-g gtf_file]      = gtf_infile                                      [required]\n" \
"   [-r]               = region, default is whole genome.(chr1:1-20000)  [option]\n" \
"   [-l]               = bool variant,get the log2 value of depth or not [option]\n" \
"   [-O OUTPUT_FILE]   = OUTPUT file. default is 'out'                   [option]\n" \
"   [-h]               = This helpful help screen.                       [option]\n" \
"   Infiles            = bam format input file(s),at least 1 bam file.   [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	printf("%s",buffer);
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
	if (c->tid < 0) return 0;
//	int i, l;
//	for (i = l = 0; i < c->n_cigar; ++i) {
//		int op = cigar[i]&0xf;
//		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
//			l += cigar[i]>>4;
//	}
//	fprintf(stderr,"%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%c\t%s\n",
//		databuf->fp->header->target_name[c->tid],c->pos, bam_calend(c,cigar) ,c->pos + l,l,bam_cigar2qlen(c,cigar), 
//		bam1_qname(b), c->qual, c->bin,c->mtid,c->mpos,c->isize,(c->flag&BAM_FREVERSE)? '-' : '+',bam_aux2Z(bam_aux_get(b,"MD")));
	unsigned short current_GC =0;
	char *qseq=cal_GC(bam1_seq(b),c->l_qseq,&current_GC);
	Region *R=(Region *)malloc(sizeof(Region));
	fill_item(b,c,cigar,qseq,R,databuf);
//	printf("%s\n",((Region *)list_tail(*(databuf->my_list+c->tid)))->qname);
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
	r->MD=strdup(bam_aux2Z(bam_aux_get(b,"MD")));
	list_push_back(*(databuf->my_list+c->tid), r);
}

void parse_my_list(DataBuf *databuf,unsigned int j){
	databuf->max_depth=0;
	databuf->max_wig=0;
	if (!list_empty(*(databuf->my_list+j))) {
		*(databuf->wig_list+j)=list_create(Wig_delete_int);
		databuf->max_wig=tity_sam2wig(*(databuf->my_list+j),*(databuf->wig_list+j));
		databuf->Stack=(struct list **)calloc(MAX_DEPTH,sizeof(struct list *));
		databuf->max_depth = Stack_list(*(databuf->my_list+j),databuf->Stack);
	}
	else{
		list_delete(*(databuf->my_list+j));
	}
}

void Wig_delete_int(void *DATA) {
	free(DATA);
}

double tity_sam2wig(struct list *LIST,struct list *wig_list){
	struct item *current, *next;
	int i=0,l=0,j=0,prevkey=0;
	double count=0;
	HASHTBL *start = hashtbl_create((HSIZE)2*LIST->size, NULL);
	HASHTBL *end = hashtbl_create((HSIZE)2*LIST->size, NULL);
	char *buf=(char *)malloc(128*sizeof(char));
	Region *r;
	current = LIST->head;
	while (current) {
		next = current->next;
		r = (Region *)current->data;
		if (current) {
			unsigned int temp_start= r->pos;
			for (j = 0,l = 0; j < r->n_cigar; ++j) {
				int op = r->cigar[j]&0xf;
				if ( op == BAM_CDEL || op == BAM_CREF_SKIP){
					l = r->cigar[j]>>4;
					temp_start+=l;
				}
				else if (op == BAM_CMATCH ) {
					l = r->cigar[j]>>4;
					memset(buf,0,128*sizeof(char));
					sprintf(buf,"%d",temp_start);
					void *data=hashtbl_get(start,buf);
					if (data==NULL) {
						int *init=(int *)malloc(sizeof(int));
						*init=1;
						hashtbl_insert(start, buf, (void *)init);
					}
					else{
						(*(int *)data)++;
					}
					
					temp_start+=l;
					memset(buf,0,128*sizeof(char));
					sprintf(buf,"%d",temp_start);

					data=hashtbl_get(end,buf);
					if (data==NULL) {
						int *init=(int *)malloc(sizeof(int));
						*init=1;
						hashtbl_insert(end, buf, (void *)init);
					}
					else{
						(*(int *)data)++;
					}
				}
			}
		}
		current = next ;
	}
	free(buf);
	char **all_keys=union_hashed_keys(start,end);
	int all_keys_count=start->count+end->count;
	double max_wig=0;
	for (i = 0; i < all_keys_count; ++i) {
		int pos=atoi(all_keys[i]);
		if (pos==prevkey) continue;
		if (prevkey) {
			Wig *wig=(Wig *)malloc(sizeof(Wig));
			wig->start=prevkey;
			wig->end=pos;
//			fprintf(stderr,"%d\t%d\t%f\n",prevkey,pos,count);
			if (globalArgs.log){
				wig->depth=count>0?(count>1?log(count)/log(globalArgs.log):0.5):0;
			}
			else{
				wig->depth=count;
			}
			max_wig=wig->depth>max_wig?wig->depth:max_wig;
			if (!list_empty(wig_list)) {
				Wig *wig_last=(Wig *)list_tail(wig_list);
				if (wig_last->end==prevkey && wig_last->depth == wig->depth) {
					wig_last->end=pos;
					free(wig);
				}
				else{
					list_push_back(wig_list, wig);
				}
			}
			else{
				list_push_back(wig_list, wig);
			}
		}
		prevkey=pos;
		void *data=hashtbl_get(start,all_keys[i]);
		if (data!=NULL) count+= *(int*)data;
		data=hashtbl_get(end,all_keys[i]);
		if (data!=NULL) count-= *(int*)data;
	}
	hashtbl_destroy(start,NULL);
	hashtbl_destroy(end,NULL);
	return max_wig;
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
	printf("i: %d\tlist_size: %d\n",i,list_size(LIST));
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
//	printf("i: %d\tlist_size: %d\nDone stack list\n",i,list_size(LIST));
	if (!list_empty(LIST)) {
		list_delete(LIST);
	}else{
		free(LIST);
	}
	return i;
}

void creat_pdf(DataBuf *databuf,unsigned short numInfiles,unsigned int *max_transcript_num,int beg,int end,int *Width,int *Height) {
	float scale=10;
	int i=0,each_transcript_height=20;
	*Width = (end-beg)/scale>SCALE_WIDTH?SCALE_WIDTH:(end-beg)/scale + 2*LEFT;
	*Height = TOP;
	(*Height)+=(*max_transcript_num)*(each_transcript_height+5)+(*max_transcript_num>0?TOP:0);
	for (i=0;i<numInfiles ;i++ ) {
		if (globalArgs.pileup)(*Height)+=1.2*EACH_READ_HEIGHT*(databuf+i)->max_depth+((databuf+i)->max_depth>0?TOP:0);
		if (globalArgs.depth) (*Height)+=1.2*EACH_WIG_HEIGHT*(databuf+i)->max_wig+((databuf+i)->max_depth>0?TOP:0);
//		printf("max_depth: %d\tmax_wig: %f\t%f\n",(databuf+i)->max_depth,(databuf+i)->max_wig,1.2*EACH_WIG_HEIGHT*(databuf+i)->max_wig+((databuf+i)->max_wig>0?TOP:0));
	}
}

void draw_GeneList_rectangle(unsigned short transcript_num,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,int beg,int end){
	float scale=10;
	unsigned int each_transcript_height=20;
	unsigned int Height = transcript_num*(each_transcript_height+5);

	cairo_set_line_width(cr, 1);
	cairo_set_source_rgb(cr,0,0,0);
	cairo_rectangle(cr, LEFT,*current_top,(end-beg)/scale>SCALE_WIDTH?SCALE_WIDTH:(end-beg)/scale,Height);
	cairo_stroke (cr);
	(*current_top)+=Height+TOP;
}

void draw_GeneList(GeneCluster *gene_cluster,unsigned short transcript_num,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,int beg,int end) {
	struct item *current, *next;

	unsigned int each_transcript_height=20;
	float scale=10;
	float ratio= (end-beg)/scale>SCALE_WIDTH ? (float)SCALE_WIDTH/(end-beg) : 1/scale;

	float x=LEFT;
	unsigned int y=*current_top;
	current = (gene_cluster->GeneList)->head;
	while (current) {
		y+=5;
		x=LEFT;
		next=current->next;
		Gene *current_gene=(Gene *)current->data;
		int *current_color=(((GtfLine *)current_gene->gtflist->head->data)->strand == '+' ? rgb[1] : rgb[10]);
//		printf("%u\t%u\n",current_gene->start,current_gene->end);
		cairo_set_line_width(cr, 0.1);
		cairo_set_source_rgb(cr,(double)current_color[0]/255,(double)current_color[1]/255,(double)current_color[2]/255);
		struct item *current2,*next2;
		current2=current_gene->gtflist->head;
		while (current2) {
			next2=current2->next;
			GtfLine *gtfline=(GtfLine *)current2->data;
			if (!strcmp(gtfline->feature,"exon")) {
				cairo_rectangle(cr, (int)(x+(gtfline->start-beg)*ratio),y,(int)(gtfline->end-gtfline->start)*ratio,each_transcript_height);
				cairo_fill(cr);
//				fprintf(stderr,"%s\t%s\t%u\t%u\t%c\t%c\t%s\n",gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
			}
			current2=next2;
		}
		cairo_move_to(cr, (int)(x+(current_gene->start-beg)*ratio),y+each_transcript_height/2);
		cairo_line_to (cr, (int)(x+(current_gene->end-beg)*ratio),y+each_transcript_height/2);
		cairo_stroke (cr);
		y+=each_transcript_height;
		current=next;
	}
}

void draw_wig(struct list *wig_list,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,double max_wig,int beg,int end,const char *filename) {
	if (!wig_list) return;
	struct item *current, *next;
	float scale=10;
	float ratio= (end-beg)/scale>SCALE_WIDTH ? (float)SCALE_WIDTH/(end-beg) : 1/scale;
	unsigned int Height = 1.2*max_wig*EACH_WIG_HEIGHT;
	cairo_set_line_width (cr, 0.1);
	cairo_set_source_rgb(cr,0,0,0);
	cairo_rectangle(cr, LEFT,*current_top, (end-beg)/scale>SCALE_WIDTH?SCALE_WIDTH:(end-beg)/scale, Height);
	cairo_stroke (cr);
	
	float x=LEFT;
	unsigned int y=*current_top;
	int *current_color=rgb[3];

	current = wig_list->head;
	Wig *wig = (Wig *)current->data;
	cairo_move_to(cr, (int)(x+(wig->start-beg+1)*ratio), (int)y);
	while (current) {
		next = current->next;
		wig = (Wig *)current->data;
		if (current) {
			cairo_line_to(cr, (int)(x+(wig->start-beg+1)*ratio),(int)y + wig->depth * EACH_WIG_HEIGHT);
			if (!wig->depth){
				cairo_close_path (cr);
				cairo_set_source_rgb(cr,(double)current_color[0]/255,(double)current_color[1]/255,(double)current_color[2]/255);
				cairo_fill(cr);
				cairo_move_to(cr,(int)(x+(wig->end-beg+1)*ratio), (int)y + wig->depth * EACH_WIG_HEIGHT);
			}
			else{
				cairo_line_to (cr, (int)(x+(wig->end-beg+1)*ratio), (int)y + wig->depth * EACH_WIG_HEIGHT);
			}
		}
		current = next ;
	}
	cairo_line_to (cr, (int)(x+(wig->end-beg+1)*ratio), (int)y);
	cairo_close_path (cr);
	cairo_set_source_rgb(cr,(double)current_color[0]/255,(double)current_color[1]/255,(double)current_color[2]/255);
	cairo_fill(cr);
				
	list_delete(wig_list);
	(*current_top)+=Height+TOP;
}

void draw_region(struct list **Stack,cairo_surface_t *surface,cairo_t *cr,unsigned int *current_top,unsigned int max_depth,int beg,int end,const char *filename) {
	unsigned int i=0,j=0,l=0;
	struct item *current, *next;
	Region *r;

	unsigned int Height = 1.2*max_depth*EACH_READ_HEIGHT;
	float scale=10;
	float ratio= (end-beg)/scale>SCALE_WIDTH ? (float)SCALE_WIDTH/(end-beg) : 1/scale;
	
	cairo_set_line_width (cr, 1);
	cairo_set_source_rgb(cr,0,0,0);
	cairo_rectangle(cr, LEFT,*current_top, (end-beg)/scale>SCALE_WIDTH?SCALE_WIDTH:(end-beg)/scale, Height);
	cairo_stroke (cr);

	float x=LEFT;
	unsigned int y=*current_top;
	cairo_set_font_size (cr, 8);
	cairo_select_font_face (cr, "Sans",CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_move_to(cr, x,y-TOP/2);
	cairo_show_text(cr, filename);
	
	const double styleDotted[2]={1,1};
	char *SNP=(char *)malloc(2*sizeof(char));
	for (i=0;i<max_depth ;i++ ) {
		x=LEFT;
		y+=EACH_READ_HEIGHT;
		if (list_empty(*(Stack+i))) break;
		current = (*(Stack+i))->head;
		while (current) {
			next = current->next;
			r = (Region *)current->data;
			if (current) {
				int *current_color= r->strand == '+' ? rgb[1] : rgb[10];
				unsigned int temp_start= r->pos;
				char *pEnd=r->MD;
				long SNP_pos=0;
				for (j = 0,l = 0; j < r->n_cigar; ++j) {
					int op  = bam_cigar_op(r->cigar[j]);
					int len = bam_cigar_oplen(r->cigar[j]);
					if ( op == BAM_CREF_SKIP){
						cairo_set_line_width (cr, 0.1);
						cairo_set_dash (cr,styleDotted,2,0);
						cairo_set_source_rgba (cr,(double)current_color[0]/255,(double)current_color[1]/255,(double)current_color[2]/255,0.4);
						cairo_move_to(cr, (int)(x+(temp_start-beg+1)*ratio), (int)y);
						cairo_line_to (cr, (int)(x+(temp_start+len-beg+1)*ratio), (int)y);
						cairo_stroke (cr);
						temp_start+=len;
						SNP_pos+=len;
						l += len;
					}
					else if (op == BAM_CDEL ) {
						cairo_set_line_width (cr, 1);
						cairo_set_dash (cr,styleDotted,0,0);
						cairo_set_source_rgba (cr,(double)rgb[4][0]/255,(double)rgb[4][1]/255,(double)rgb[4][2]/255,1);
						cairo_move_to(cr, (int)(x+(temp_start-beg+1)*ratio), (int)y);
						cairo_line_to (cr, (int)(x+(temp_start+len-beg+1)*ratio), (int)y);
						cairo_stroke (cr);
						
						cairo_set_font_size (cr, 20);
						cairo_set_source_rgb (cr, 0, 0, 0);
						cairo_move_to(cr, (int)(x+(temp_start-beg+1)*ratio), (int)y);
						SNP[0]=bam_cigar_opchr(r->cigar[j]);SNP[1]=0;
						cairo_show_text(cr,SNP );
						temp_start+=len;
						SNP_pos+=len;
						l += len;
					}
					else if (op == BAM_CMATCH) {
						cairo_set_line_width (cr, 2);
						cairo_set_dash (cr,styleDotted,0,0);
						cairo_set_source_rgb(cr,(double)current_color[0]/255,(double)current_color[1]/255,(double)current_color[2]/255);
						cairo_move_to(cr, (int)(x+(temp_start-beg+1)*ratio), (int)y);
						cairo_line_to (cr, (int)(x+(temp_start+len-beg+1)*ratio), (int)y);
						cairo_stroke (cr);
						while (SNP_pos>=l && SNP_pos<l+len) {
							SNP_pos+=strtol(pEnd,&pEnd,10);
							if (! *pEnd || SNP_pos>l+len) break;
//							if (SNP_pos>100)
//								printf("%s\t%d\t%d\t%d\t%s\t%s\t%c\t%s\t%ld\t%u\t%d\n",r->chromosome,r->pos,r->qlen,r->end,r->qname,r->qseq,r->strand,r->MD,SNP_pos,l,len);
							SNP[0]=*pEnd;SNP[1]=0;
							while (*pEnd < '0' || *pEnd > '9'){
								if (*pEnd!='^'){
									SNP_pos++;
								}else{
									while (*pEnd < '0' || *pEnd > '9') pEnd++;
								}
								pEnd++;
							}
							cairo_set_font_size (cr, 1);
							cairo_set_source_rgb (cr, 0, 0, 0);
							cairo_move_to(cr, (int)(x+(r->pos+SNP_pos-beg+1)*ratio), (int)y);
							cairo_show_text(cr, SNP);
							
							cairo_set_source_rgb(cr,(double)rgb[5][0]/255,(double)rgb[5][1]/255,(double)rgb[5][2]/255);
							cairo_set_line_width (cr, 0.1);
							cairo_move_to(cr, (int)(x+(r->pos+SNP_pos-beg+1)*ratio), (int)y-EACH_READ_HEIGHT/4);
							cairo_line_to (cr, (int)(x+(r->pos+SNP_pos-beg+1)*ratio), (int)y+EACH_READ_HEIGHT/4);
							cairo_stroke (cr);
						}
						temp_start+=len;
						l += len;
					}
				}
//				cairo_set_font_size (cr, 1);
//				cairo_set_source_rgb (cr, 0, 0, 0);
//				cairo_move_to(cr, (int)(x+(r->pos-beg+1)*ratio), (int)y+EACH_READ_HEIGHT/4);
//				cairo_show_text(cr, r->qseq);
//				printf("%d\t%d\t%d\t%d\t#%d#\n",i,r->pos,r->end,(*(Stack+i))->size,r->n_cigar);
				item_delete(*(Stack+i), current);
				(*(Stack+i))->size--;
			}
			current = next ;
		}
		free(*(Stack+i));
	}
	(*current_top)+=Height+TOP;
	free(Stack);
	free(SNP);
}

char *create_pdf(const char *chromosome,int beg,int end,const char *filename){
	char *output=(char *)malloc(1024*sizeof(char));
	memset(output,0,1024*sizeof(char));
	sprintf(output,"%s_%s_%d_%d.pdf",filename,chromosome,beg,end);
	return output;
}

void free_cario(cairo_t *cr,cairo_surface_t *surface,char *filename){
	cairo_destroy (cr);
	cairo_surface_destroy (surface);
	free(filename);
}

void Region_delete_int(void *DATA) {
	Region *R=(Region *)DATA;
	free(R->cigar);
	free(R->chromosome);
	free(R->qname);
	free(R->qseq);
	free(R->MD);
	free(R);
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

int test_bsearch(HASHTBL *Clusterd_GTF,DataBuf *databuf,int ref,int beg,int end) {
	ChrCluster *chr_culster=(ChrCluster *)hashtbl_get(Clusterd_GTF, databuf->fp->header->target_name[ref]);
	int gene_cluster_index=bsearch_04(chr_culster,beg);
	if (gene_cluster_index>0) {
		printf("gene_cluster_index: %d\tgene_cluster_start: %u\tgene_cluster_end: %u\n",gene_cluster_index,chr_culster->chr_gene_cluster[gene_cluster_index].start,chr_culster->chr_gene_cluster[gene_cluster_index].end);
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
		printf("max_transcript_num: %u\n",*max_transcript_num);
	}
	return gene_cluster_index;
}

void draw_a_region(HASHTBL *Clusterd_GTF,DataBuf *databuf){
	int ref, beg, end,i,boolen=0;
	unsigned int current_top=50;
//	printf("%s\t%d\t%d\n",databuf->fp->header->target_name[ref],beg,end);
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		ref=0;beg=0;end=0;
		bam_parse_region((databuf+i)->fp->header, globalArgs.region, &ref, &beg, &end);
		if (ref < 0) {
			printf("bam2bed: Invalid region %s\n", globalArgs.region);
			exit(1);
		}
		*((databuf+i)->my_list+ref)=list_create(Region_delete_int);
		bam_fetch((databuf+i)->fp->x.bam, (databuf+i)->idx, ref, beg, end, (databuf+i), fetch_func);
		printf("%s read count in this region: %d\n",*(globalArgs.infiles+i),list_size(*((databuf+i)->my_list+ref)));
		if (!list_empty(*((databuf+i)->my_list+ref))) boolen=1;
		parse_my_list(databuf+i,ref);
	}
	if (boolen) {
		unsigned int max_transcript_num=0;
//		int gene_cluster_index = test_bsearch(Clusterd_GTF,databuf,ref,beg,end);
//		printf("Done binary search GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
		ChrCluster *chr_culster=(ChrCluster *)hashtbl_get(Clusterd_GTF, databuf->fp->header->target_name[ref]);
		int gene_cluster_index=bsearch_draw(chr_culster,beg,end,&max_transcript_num);
		
		char *outfile=create_pdf(databuf->fp->header->target_name[ref],++beg,end,globalArgs.outfile);
		int Width=0;
		int Height=0;
		creat_pdf(databuf,globalArgs.numInfiles,&max_transcript_num, beg,end,&Width, &Height);
		cairo_surface_t *surface = cairo_pdf_surface_create(outfile, Width, Height);
		cairo_t *cr = cairo_create(surface);
		if (gene_cluster_index>0) {
			draw_GeneList(&chr_culster->chr_gene_cluster[gene_cluster_index],max_transcript_num,surface,cr,&current_top,beg,end);
			while (chr_culster->chr_gene_cluster[gene_cluster_index].end <end) {
				draw_GeneList(&(chr_culster->chr_gene_cluster[++gene_cluster_index]),max_transcript_num,surface,cr,&current_top,beg,end);
			}
			draw_GeneList_rectangle(max_transcript_num,surface,cr,&current_top,beg,end);
		}
		for (i=0;i<globalArgs.numInfiles ;i++ ) {
			if ((databuf+i)->max_depth>0) {
				if (globalArgs.pileup) draw_region((databuf+i)->Stack,surface,cr,&current_top,(databuf+i)->max_depth,beg,end,*(globalArgs.infiles+i));
				if (globalArgs.depth) draw_wig(*((databuf+i)->wig_list+ref),surface,cr,&current_top,(databuf+i)->max_wig,beg,end,*(globalArgs.infiles+i));
			}
		}
		free_cario(cr,surface,outfile);
	}
}

void draw_all_genes(HASHTBL *Clusterd_GTF,DataBuf *databuf,const char *current_dir,char *work_dir){
	unsigned int i,j,m=0;
	for (j=0;j<databuf->fp->header->n_targets ;j++ ) {
		sprintf(work_dir,"%s/%s",current_dir,databuf->fp->header->target_name[j]);
		mkdir(work_dir,0750);
		int status=chdir(work_dir);
		if (status) exit(0);
		ChrCluster *chr_culster=(ChrCluster *)hashtbl_get(Clusterd_GTF, databuf->fp->header->target_name[j]);
		
		if (chr_culster==NULL) continue;
		printf("%s\t%u\n",databuf->fp->header->target_name[j],chr_culster->Cluster_count);
		for (m=0;m<chr_culster->Cluster_count ;m++ ) {
			int boolen=0;
			for (i=0;i<globalArgs.numInfiles ;i++ ) {
				*((databuf+i)->my_list+j)=list_create(Region_delete_int);
				fprintf(stderr,"%d\t%d\n",chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end);
				bam_fetch((databuf+i)->fp->x.bam, (databuf+i)->idx,j,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,(databuf+i), fetch_func);
				printf("%u\t%s read count in this region: %d\n",m+1,*(globalArgs.infiles+i),list_size(*((databuf+i)->my_list+j)));
				if (!list_empty(*((databuf+i)->my_list+j))) boolen=1;
				parse_my_list(databuf+i,j);
			}
			if (boolen) {
				unsigned int max_transcript_num=list_size(chr_culster->chr_gene_cluster[m].GeneList);
				unsigned int current_top=50;

				int Width=0;
				int Height=0;
				char *outfile=create_pdf(databuf->fp->header->target_name[j],chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,globalArgs.outfile);
				creat_pdf(databuf,globalArgs.numInfiles,&max_transcript_num, chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end, &Width, &Height);
				cairo_surface_t *surface = cairo_pdf_surface_create(outfile, Width, Height);
				cairo_t *cr = cairo_create(surface);

				draw_GeneList(&chr_culster->chr_gene_cluster[m],max_transcript_num,surface,cr,&current_top,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end);
				draw_GeneList_rectangle(max_transcript_num,surface,cr,&current_top,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end);
				
				for (i=0;i<globalArgs.numInfiles ;i++ ) {
					if ((databuf+i)->max_depth>0) {
						if (globalArgs.pileup) draw_region((databuf+i)->Stack,surface,cr,&current_top,(databuf+i)->max_depth,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,*(globalArgs.infiles+i));
						if (globalArgs.depth) draw_wig(*((databuf+i)->wig_list+j),surface,cr,&current_top,(databuf+i)->max_wig,chr_culster->chr_gene_cluster[m].start,chr_culster->chr_gene_cluster[m].end,*(globalArgs.infiles+i));
					}
				}
				free_cario(cr,surface,outfile);
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
	globalArgs.depth=0;
	globalArgs.pileup=0;
	globalArgs.region="-";
	globalArgs.outfile="out";
	globalArgs.log=2;
	globalArgs.gtf="-";

	const char *optString = "o:DPg:r:l:O:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.overlap = atoi(optarg);
				break;
			case 'D':
				globalArgs.depth = 1;
				break;
			case 'P':
				globalArgs.pileup = 1;
				break;
			case 'r':
				globalArgs.region = optarg;
				break;
			case 'l':
				globalArgs.log = atoi(optarg);
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
				printf("error parameter!\n");
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	globalArgs.infiles = argv + optind;
	globalArgs.numInfiles = argc - optind;

	long long begin;
	begin=usec();

	HASHTBL *GTF=load_gtf(globalArgs.gtf);
	printf("Done load GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
//	print_gtf(GTF);
	HASHTBL *Clusterd_GTF=cluster_gtf(GTF);
	printf("Done clusted GTF at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));

	unsigned int i=0;
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			printf( "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
			exit(1);
		}
		(databuf+i)->my_list = (struct list **)calloc((databuf+i)->fp->header->n_targets,sizeof(struct list *));
		(databuf+i)->wig_list = (struct list **)calloc((databuf+i)->fp->header->n_targets,sizeof(struct list *));
		if (((databuf+i)->idx = bam_index_load(*(globalArgs.infiles+i))) == 0) {
			printf( "bam2bed: BAM indexing file is not available.\n");
			exit(1);
		}
//		printf("$%d\t%d\n",j,(databuf+i)->fp->header->target_len[j]);
	}

	char *current_dir=(char *)calloc(128,sizeof(char));
	current_dir = getcwd(current_dir, 128*sizeof(char));
	char *work_dir=(char *)calloc(256,sizeof(char));
	
	if (strncmp(globalArgs.region,"-",6)==0) { /* if a region is not specified */
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
		free((databuf+i)->wig_list);
		samclose((databuf+i)->fp);
	}
	free(databuf);
	free(current_dir);
	free(work_dir);

	printf("Finished at %.3f sec\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
