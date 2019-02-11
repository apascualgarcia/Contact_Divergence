/* 
   Program Contact_divergence
   Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)
   ubastolla@cbm.uam.es
   Reads a multiple alignment and computes contact divergence and other
   structure comparison measures.

   INPUT: file with multiple alignment in FASTA format (at least 2)
   Protein names must be names of PDB files.
   The first line may be PDBDIR=<directory of PDB files>
   (default: current directory)

   OUTPUT: For each protein pair, structural scores printed are
   Contact_divergence, contact overlap, TM score (default no)

*/

#include "Contact_divergence_aux.h"
#include "D_Cont.h"
#include "protein.h"
#include "cont_list.h"
#include "allocate.h"
#include "read_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EXT_OUT ".sim"  // Extension for output file

#define IJ_MIN_DEF 3
#define CONT_TYPE_DEF 'c'
#define CONT_THR_DEF 4.5
int IJ_MIN=IJ_MIN_DEF;        // Only contacts with |i-j|>=IJ_MIN
float CONT_THR=CONT_THR_DEF;
char CONT_TYPE=CONT_TYPE_DEF;  //a=alpha b=beta c=all atoms
char CONT_STRING[80];
char CODENAME[40]="main_Contact_divergence.c";
int CONT_DEF=1;

struct Prot_input{
  char name[80];
  char chain;
  char *seq;
};

void help(char *pname);
void Get_input(char *file_ali, int *PRINT_TM, int *PRINT_CONT, int *PRINT_SEQ,
	      int argc, char **argv);
int Get_alignment(char ***Prot_name, char **Prot_chain, char ***Prot_seq,
		  int *Nali, char *PDB_PATH, char *file_ali);
void Set_contact_type();


/*********************************************************************
                          MAIN routine
**********************************************************************/
int main(int argc, char **argv)
{
  // INPUT
  printf("Starting %s\n", argv[0]);
  char file_ali[200];
  int PRINT_TM=0, PRINT_CONT=0, PRINT_SEQ=0;
  Get_input(file_ali, &PRINT_TM, &PRINT_CONT, &PRINT_SEQ, argc, argv);
  Set_contact_type();

  char **Prot_name, **Prot_seq, *Prot_chain, PDB_PATH[200]="./"; //
  int N_ali;
  int N_prot=Get_alignment(&Prot_name, &Prot_chain, &Prot_seq, &N_ali,
			   PDB_PATH, file_ali);

  /**************+++   READ PROTEINS  ********************/
  // Look for contact matrix and sequence files. If not found,
  // read coordinates from PDB files and print contact matrices
  struct protein *prots=malloc(N_prot*sizeof(struct protein)), *prot=prots;
  int N_pdb=0, i, *index=malloc(N_prot*sizeof(int));
  int **Prot_ali=Allocate_mat2_i(N_prot, N_ali);
    printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	   CONT_TYPE, CONT_THR,IJ_MIN);
  for(i=0; i<N_prot; i++){
    if(Read_PDB_compress(prot, Prot_name[i], &(Prot_chain[i]), PDB_PATH)>0){
      if(Align_seq(Prot_ali[N_pdb],N_ali,Prot_seq[i],prot->aseq,prot->len)<0)
	continue;
      int NC=Compute_contact_list(prot, CONT_TYPE, CONT_THR, IJ_MIN);
      printf("%d contacts\n", NC);
      index[N_pdb]=i; prot++; N_pdb++;
    }
  }
  printf("%d proteins read out of %d listed in %s\n", N_pdb, N_prot, file_ali);
  if(N_pdb<2){
    printf("ERROR, zero protein pairs found\n"); exit(8);
  } 

  
  // Pairwise computations
  char FILE_OUT[100];
  Change_ext(FILE_OUT, file_ali, EXT_OUT);
  FILE *file_out=fopen(FILE_OUT, "w");
  fprintf(file_out, "# Protein pair Cont_Divergence relatedness");
  if(PRINT_CONT)fprintf(file_out, " Cont_Overlap");
  if(PRINT_SEQ)fprintf(file_out, " Seq_Id");
  fprintf(file_out, "\n");
  for(i=0; i<N_pdb; i++){
    struct protein *proti=prots+i, *protj=prots; int j, homo;
    for(j=0; j<i; j++){
      float Seqid=Seq_identity(Prot_seq[index[i]], Prot_seq[index[j]], N_ali);
      float overlap=Compute_overlap(proti->Cont_map, proti->len, Prot_ali[i],
				    protj->Cont_map, protj->len, Prot_ali[j],
				    N_ali);
      float D_Cont=Compute_Dcont(overlap, proti->len, protj->len, &homo);
      fprintf(file_out, "%s\t%s", proti->name, protj->name);
      fprintf(file_out, "\t%.3f\t%d", D_Cont, homo);
      if(PRINT_CONT)fprintf(file_out, "\t%.3f", overlap);
      if(PRINT_SEQ)fprintf(file_out, "\t%.3f", Seqid);
      /*if(PRINT_TM){
	float tm=TM_SCORE(proti->xca, proti->len, Prot_ali[i],
	protj->xca, protj->len, Prot_ali[j], N_ali);
	}*/
      fprintf(file_out, "\n");
      protj++;
    }
  }
  printf("Results written in file %s\n", FILE_OUT);
  return(0);
}

void Get_input(char *file_ali, int *PRINT_TM, int *PRINT_CONT, int *PRINT_SEQ,
	       int argc, char **argv)
{
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);
  strcpy(file_ali, argv[1]);
  printf("Alignment file: %s\n", file_ali);
  int i;
  for(i=2; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){
      help(argv[0]);
      //}else if(strncmp(argv[i], "-tm", 3)==0){
      //*PRINT_TM=1;
    }else if(strncmp(argv[i], "-cont", 5)==0){
      *PRINT_CONT=1;
    }else if(strncmp(argv[i], "-seq", 5)==0){
      *PRINT_SEQ=1;
    }else{
      printf("WARNING, unrecognized option %s\n", argv[i]);
    }
  }
}

void help(char *pname){
  printf("Program %s\n", pname);
  printf("Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa\n");
  printf("Email: <ubastolla@cbm.uam.es>\n\n");
  printf("INPUT: Multiple sequence alignment of proteins in FASTA format\n");
  printf("The protein name is the name of a PDB file, optionally followed\n");
  printf("by the chain index (Ex: >1opd.pdb A)\n");
  printf("The first line may be PDBDIR=<directory of PDB files>\n");
  printf("(default: current directory)\n\n");
  printf("OUTPUT: For each pair proteins, structural scores are printed:\n");
  printf("Contact divergence (Pascual-García et al Proteins 2010 78:181-96)\n");
  printf("Contact overlap (optional)\n");
  printf("TM score (Zhang & Skolnick Proteins 2004 57:702-10, optional)\n\n");
  printf("OPTIONS:\n");
  //printf("-tm Print TM score\n");
  printf("-cont Print contact overlap\n");
  printf("-seq  Print sequence identity\n");
  printf("\n");
  exit(8);
}

int Get_alignment(char ***Prot_name, char **Prot_chain, char ***Prot_seq,
		  int *Nali, char *PDB_PATH, char *file_ali)
{
  // Open file
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int dir=0, n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if((n==0)&&(strncmp(string, "PDBDIR", 6)==0)){
      sscanf(string+7,"%s", PDB_PATH);
      printf("Directory for PDB files: %s\n", PDB_PATH); dir=1;
    }else if(string[0]=='>'){
      n++;
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  // Allocate and read sequences
  int CNAME=50, LMAX=10000, l=0, i;
  char chain[10], dumm[40];
  char *Seq=malloc(LMAX*sizeof(char)), *s=NULL;
  *Prot_seq=malloc(n*sizeof(char *));
  *Prot_name=malloc(n*sizeof(char *));
  *Prot_chain=malloc(n*sizeof(char));
  n=0; *Nali=0;
  file_in=fopen(file_ali, "r");
  if(dir)fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if((n==0)&&(strncmp(string, "PDBDIR", 6)==0))continue;
    if(string[0]=='>'){
      (*Prot_name)[n]=malloc(CNAME*sizeof(char));
      sscanf(string+1, "%s", (*Prot_name)[n]);
      for(i=0; i<10; i++)chain[i]='\0';
      int c=sscanf(string, "%s%s\n", dumm, chain);
      if((c>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){(*Prot_chain)[n]=chain[0];}
      else{(*Prot_chain)[n]=' ';}
      printf("%s %c\n", (*Prot_name)[n], (*Prot_chain)[n]);
      if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_seq)[0]=malloc(*Nali*sizeof(char));
	s=(*Prot_seq)[0];
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_seq)[n]=malloc(*Nali*sizeof(char));
	s=(*Prot_seq)[n]; l=0;
      }else{
	s=Seq; l=0;
      }
      n++;
    }else{
      char *c=string;
      while(*c!='\n'){*s=*c; l++; s++; c++;}
      if(l > LMAX){
	printf("ERROR, alignment length larger than maximum allowed %d\n", l);
	printf("Increase LMAX in code %s\n", CODENAME); exit(8);
      }
      if((*Nali)&&(l>*Nali)){
	printf("ERROR, too many column in alignment %d.",n);
	printf(" Expected %d, found >= %d\n", *Nali, l); exit(8); 
      }
    }
  }
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
}

void Set_contact_type(){

  if(CONT_TYPE=='a'){strcpy(CONT_STRING, "Alpha");}
  else if(CONT_TYPE=='b'){strcpy(CONT_STRING, "Beta");}
  else if(CONT_TYPE=='c'){strcpy(CONT_STRING, "All atoms");}
  else{
    printf("WARNING, undefined contact %c\n", CONT_TYPE);
    CONT_TYPE='c'; strcpy(CONT_STRING, "All atoms");
    printf("Using default %s\n", CONT_STRING);
  }
  // Default type of contacts?
  if((CONT_TYPE!=CONT_TYPE_DEF)||(CONT_THR!=CONT_THR_DEF)||
     (IJ_MIN!=IJ_MIN_DEF))CONT_DEF=0;
}
