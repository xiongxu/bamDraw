#include <stdint.h>
#include <string.h>
#include "hashtbl.h"
#include "hash-master/list.h"
#define MAX_CHR 64
#define MAX_GENE_NUM 30000

typedef struct _gtfLine_
{
	char *seqname;
	char *feature;
	unsigned int start;
	unsigned int end ;
	char strand;
	char frame;
	char *attributes;
} GtfLine;

typedef struct _Gene_
{
	unsigned int start;
	unsigned int end ;
	struct list *gtflist;
} Gene;

typedef struct _Chr_
{
	unsigned int gene_count;
	Gene *chr_gene;
} Chr;

typedef struct _GeneCluster_
{
	unsigned int start;
	unsigned int end ;
	struct list *GeneList;
} GeneCluster;

typedef struct _ChrCluster_
{
	unsigned int Cluster_count;
	GeneCluster *chr_gene_cluster;
} ChrCluster;

/*load GTF to a hash table*/
static inline void data_delete_listT(void *DATA) ;
static inline GtfLine *readNextGtf(FILE *f,const char *buf);
HASHTBL *load_gtf(const char *infile);
static inline void assign_gene_start_end(Gene *gene);
static inline int cmp_gene(const void *a, const void *b);
static inline void print_GtfList(struct list *gtfline_list);
static inline void print_gtf(HASHTBL *hash_gtf);

/*cluster GTF*/
static inline void data_delete_list_GeneT(void *DATA);
static inline int overlap(unsigned int S1,unsigned int E1,unsigned int S2,unsigned int E2);
static inline HASHTBL *cluster_gtf(HASHTBL *hash_gtf);
static inline void print_GeneList(struct list *GeneList);
static inline int cmp_GeneCluster(const void *a, const void *b);
void print_clusterd_gtf(HASHTBL *hash_clusterd_gtf);

void data_delete_listT(void *DATA) {
	GtfLine *gtfline=(GtfLine *)DATA;
	free(gtfline->seqname);
	free(gtfline->feature);
	free(gtfline->attributes);
	free(gtfline);
}

GtfLine *readNextGtf(FILE *f,const char *buf){
	GtfLine *gtfline=(GtfLine *)malloc(sizeof(GtfLine));
	gtfline->seqname=(char *)calloc(128,sizeof(char));
	gtfline->feature=(char *)calloc(128,sizeof(char));
	gtfline->attributes=(char *)calloc(2048,sizeof(char));
	sscanf(buf, "%s\t%*[^\t]\t%s\t%u\t%u\t%*[^\t]\t%c\t%c\t%*[^;];%*s \"%[^\"]",
		gtfline->seqname,gtfline->feature, &gtfline->start, &gtfline->end,&gtfline->strand,&gtfline->frame,gtfline->attributes);
	gtfline->seqname=realloc(gtfline->seqname,(strlen(gtfline->seqname)+1)*sizeof(char));
	gtfline->feature=realloc(gtfline->feature,(strlen(gtfline->feature)+1)*sizeof(char));
	gtfline->attributes=realloc(gtfline->attributes,(strlen(gtfline->attributes)+1)*sizeof(char));
//	fprintf(stderr,"%s\t%s\t%u\t%u\t%c\t%c\t%s\n",
//		gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
	return gtfline;
}

HASHTBL *load_gtf(const char *infile){
	FILE *f=fopen(infile,"r");
	if(f == NULL) {
		printf("Failed to open file (%s)", infile);
		exit(1);
	}
	char *buf=(char *)calloc(2048,sizeof(char));
	char *temp_seqname = "-";
	char *temp_attributes = "-";
//	struct hash *hash_gtf = hash_create(data_delete_hashT, MAX_CHR);
	HASHTBL *hash_gtf=hashtbl_create(MAX_CHR, NULL);
	unsigned int i=0;
	Gene *temp_chr_gene=(Gene *)calloc(MAX_GENE_NUM,sizeof(Gene));
	while (1) {
		buf=fgets(buf,2048*sizeof(char),f);
		if (feof(f)) break;
		GtfLine *gtfline=readNextGtf(f,buf);
		if (strcmp(temp_attributes,gtfline->attributes)) {
			if (strcmp(temp_attributes,"-")) {
				assign_gene_start_end(&temp_chr_gene[i++]);
			}
			if (strcmp(temp_seqname,gtfline->seqname)) {
//				fprintf(stderr,"%u\n",i);
				if (strcmp(temp_seqname,"-")) {
					Chr *CHR=(Chr *)malloc(sizeof(Chr));
					CHR->gene_count=i;
					CHR->chr_gene=(Gene *)calloc(CHR->gene_count,sizeof(Gene));
					memmove(CHR->chr_gene,temp_chr_gene,(CHR->gene_count)*sizeof(Gene));
//					hash_insert(hash_gtf,temp_seqname,CHR);
					hashtbl_insert(hash_gtf,temp_seqname,(void *)CHR);
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
//	hash_insert(hash_gtf,temp_seqname,CHR);
	hashtbl_insert(hash_gtf,temp_seqname,(void *)CHR);
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
		printf("%s\t%s\t%u\t%u\t%c\t%c\t%s\n",gtfline->seqname,gtfline->feature, gtfline->start, gtfline->end,gtfline->strand,gtfline->frame,gtfline->attributes);
		current = next;
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

void print_gtf(HASHTBL *hash_gtf){
	ENTRY **nodes=dump_hash_table(hash_gtf);
	qsort(nodes, hash_gtf->count, sizeof(ENTRY *), compare_hashed_data_key);

	unsigned i=0,j=0;
	fprintf(stdout,"%d\n",(int)hash_gtf->count);
	for (i=0;i<hash_gtf->count ;i++ ) {
		Chr *chr=(Chr *)nodes[i]->data;
		printf("%s: %u\n",(char *)nodes[i]->key,chr->gene_count);
		qsort(chr->chr_gene,chr->gene_count,sizeof(Gene),cmp_gene);
		for (j=0;j< chr->gene_count;j++ ) {
			Gene current_gene=chr->chr_gene[j];
			printf("%u\t%u\n",current_gene.start,current_gene.end);
			print_GtfList(current_gene.gtflist);
		}
	}
}

void data_delete_list_GeneT(void *DATA) {
	Gene *gene=(Gene *)DATA;
	list_delete(gene->gtflist);
	free(gene);
}

int overlap(unsigned int S1,unsigned int E1,unsigned int S2,unsigned int E2){
	return E1>=S2?1:0;
}

HASHTBL *cluster_gtf(HASHTBL *hash_gtf) {
//	struct hash *hash_cluster_gtf=hash_create(data_delete_hash_clusterT, MAX_CHR);
	HASHTBL *hash_cluster_gtf=hashtbl_create(MAX_CHR, NULL);
	unsigned i=0,j=0;
	ENTRY *current;
	for (i=0;i<hash_gtf->size ;i++ ) {
		current = hash_gtf->nodes[i];
		while (current) {
			Chr *chr=(Chr *)current->data;
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
//			fprintf(stderr,"%s\t%u\n",current->key,chr_culster->Cluster_count);
			chr_culster->chr_gene_cluster=(GeneCluster *)calloc(chr_culster->Cluster_count,sizeof(GeneCluster));
			memmove(chr_culster->chr_gene_cluster,chr_gene_cluster,chr_culster->Cluster_count*sizeof(GeneCluster));
//			hash_insert(hash_cluster_gtf,current->key,chr_culster);
			hashtbl_insert(hash_cluster_gtf,current->key,(void *)chr_culster);
			free(chr_gene_cluster);
			current = current->next;
		}
	}
	return hash_cluster_gtf;
}

void print_GeneList(struct list *GeneList){
	struct item *current, *next;
	if (!GeneList) return;
	printf("GeneListSize: %d\n",list_size(GeneList));
	current = GeneList->head;
	while (current) {
		next = current->next;
		Gene *current_gene=(Gene *)current->data;
		printf("Gene_start: %u\tGene_end: %u\n",current_gene->start,current_gene->end);
		print_GtfList(current_gene->gtflist);
		current = next;
	}
}

int cmp_GeneCluster(const void *a, const void *b){
	GeneCluster *c=(GeneCluster *)a;
	GeneCluster *d=(GeneCluster *)b;
	if (c->start!=d->start) {
		return c->start-d->start;
	}
	else{
		return c->end-d->end;
	}
}

void print_clusterd_gtf(HASHTBL *hash_clusterd_gtf){
	ENTRY **nodes=dump_hash_table(hash_clusterd_gtf);
	qsort(nodes, hash_clusterd_gtf->count, sizeof(ENTRY *), compare_hashed_data_key);
	unsigned i=0,j=0;
	fprintf(stdout,"%d\n",(int)hash_clusterd_gtf->count);
	for (i=0;i<hash_clusterd_gtf->count ;i++ ) {
		ChrCluster *chr_cluster=(ChrCluster *)nodes[i]->data;
		printf("%s: %u\n",(char *)nodes[i]->key,chr_cluster->Cluster_count);
		qsort(chr_cluster->chr_gene_cluster,chr_cluster->Cluster_count,sizeof(GeneCluster),cmp_GeneCluster);
		for (j=0;j< chr_cluster->Cluster_count;j++ ) {
			GeneCluster current_gene_cluster=chr_cluster->chr_gene_cluster[j];
			printf("%u\t%u\n",current_gene_cluster.start,current_gene_cluster.end);
			print_GeneList(current_gene_cluster.GeneList);
		}
	}
}

