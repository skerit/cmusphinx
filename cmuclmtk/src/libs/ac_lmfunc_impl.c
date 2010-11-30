/* ====================================================================
 * Copyright (c) 1999-2006 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

#include "ac_lmfunc_impl.h"
#include "ac_hash.h"
#include "ac_parsetext.h"
#include "idngram2lm.h"  // in liblmest/

/* Start Implementation of text2wfreq */

int text2wfreq_impl(FILE* infp, FILE* outfp, int init_nwords, int verbosity)
{
  int hash_size, scanrc;
  struct hash_table vocab;
  char word[MAX_STRING_LENGTH];

  hash_size = nearest_prime( init_nwords );
  new_hashtable( &vocab, hash_size );
  while( (scanrc = fscanf(infp, "%500s", word )) == 1 ) {
    if ( strlen( word ) >= MAX_STRING_LENGTH ) {
      pc_message(verbosity,1,"text2wfreq : WARNING: word too long, will be split: %s...\n",word);
    }
    if (strlen(word)) {
      update( &vocab, word ,verbosity);
    }
  }
  if ( scanrc != EOF ) {
    quit(-1,"Error reading input\n");
  }

  print( outfp, &vocab );
  return 0;
}

/* Start Implementation of wfreq2vocab */

typedef struct {
  char *word;
  int count;
} word_rec;

static int sort_by_count(const void *rec1,const void *rec2) {

  word_rec *r1;
  word_rec *r2;
  
  r1 = (word_rec *) rec1;
  r2 = (word_rec *) rec2;

  return(r2->count-r1->count);

}		  

static int sort_alpha(const void *rec1,const void *rec2) {
 
  word_rec *r1;
  word_rec *r2;
 
  char *s1;
  char *s2;
  
  r1 = (word_rec *) rec1;
  r2 = (word_rec *) rec2;
  
  s1 = r1->word;
  s2 = r2->word;

  return (strcmp(s1,s2));
 
}

/* To make this function less dependent on input stream, just pull records out and create an interface for it
 */
int wfreq2vocab_impl(FILE* ifp, FILE* ofp, int cutoff, int vocab_size, int num_recs, int verbosity)
{
  flag gt_set;
  flag top_set;
  int current_rec;
  int num_above_threshold;
  int num_to_output;
  int i;
  word_rec *records;
  char temp_word[750];

  gt_set = (cutoff != -1);
  top_set = (vocab_size != -1);
  if(cutoff==-1) cutoff=0;
  if(vocab_size==-1) vocab_size=0;

  if (gt_set && top_set) 
    quit(-1,"wfreq2vocab : Error : Can't use both the -top and the -gt options.\n");

  if (!gt_set && !top_set) 
    vocab_size = 20000;

  if (gt_set) 
    pc_message(verbosity,2,"wfreq2vocab : Will generate a vocabulary containing all words which\n              occurred more that %d times. Reading wfreq stream from stdin...\n",cutoff);
  else 
    pc_message(verbosity,2,"wfreq2vocab : Will generate a vocabulary containing the most\n              frequent %d words. Reading wfreq stream from stdin...\n",vocab_size);

  current_rec = 0;
  num_above_threshold = 0;

  records = (word_rec *) rr_malloc(sizeof(word_rec)*num_recs);

  while (!rr_feof(ifp)) {
    if (fscanf(ifp, "%s %d",temp_word,&(records[current_rec].count)) != 2) {
      if (!rr_feof(ifp)) 
	quit(-1,"Error reading unigram counts from standard input.\n");

    }else {
      records[current_rec].word = salloc(temp_word);
      if (gt_set && records[current_rec].count > cutoff) 
	num_above_threshold++;

      current_rec++;
    }

    if(current_rec > num_recs ){
      quit2(-1,"The number of records %d reach the user-defined limit %d, consider to increase the number of records by -records\n",current_rec,num_recs);
    }
  }

  /* Sort records in descending order of count */

  qsort((void*) records,(size_t) current_rec, sizeof(word_rec),sort_by_count);

  if (gt_set) 
    num_to_output = num_above_threshold;
  else 
    num_to_output = vocab_size;

  if (current_rec<num_to_output) 
    num_to_output = current_rec;

  /* Now sort the relevant records alphabetically */

  qsort((void*) records,(size_t) num_to_output, sizeof(word_rec),sort_alpha);

  if (gt_set) 
    pc_message(verbosity,2,"Size of vocabulary = %d\n",num_to_output);

  if (num_to_output>MAX_UNIGRAM) {
    pc_message(verbosity,1,"Warning : Vocab size exceeds %d. This might cause problems with \n",MAX_UNIGRAM);
    pc_message(verbosity,1,"other tools, since word id's are stored in 2 bytes.\n");
  }

  if (num_to_output == 0) 
    pc_message(verbosity,1,"Warning : Vocab size = 0.\n");
  /* Print the vocab to stdout */
  
  printf("## Vocab generated by v2 of the CMU-Cambridge Statistcal\n");
  printf("## Language Modeling toolkit.\n");
  printf("##\n");
  printf("## Includes %d words ",num_to_output);
  printf("##\n");

  for (i=0;i<=num_to_output-1;i++) 
    fprintf(ofp,"%s\n",records[i].word);

  pc_message(verbosity,0,"wfreq2vocab : Done.\n");

  return 0;
}

int read_vocab(char* vocab_filename, 
	       int verbosity,
	       struct idngram_hash_table* vocabulary,
	       int M
	       )
{
  FILE *vocab_file;
  int vocab_size;
  char temp_word[MAX_WORD_LENGTH];
  char temp_word2[MAX_WORD_LENGTH];

  vocab_size = 0;
  vocab_file = rr_iopen(vocab_filename);

  pc_message(verbosity,2,"Reading vocabulary... \n");

  while (fgets (temp_word, sizeof(temp_word),vocab_file)) {
    if (strncmp(temp_word,"##",2)==0) continue;
    sscanf (temp_word, "%s ",temp_word2);

    /*    printf("hey hey %s %d\n ", temp_word2, idngram_hash(temp_word2,M));*/

    /* Check for repeated words in the vocabulary */    
    if (index2(vocabulary,temp_word2) != 0)
      warn_on_repeated_words(temp_word2);

    warn_on_wrong_vocab_comments(temp_word);
    vocab_size++;
    /*    printf("%s %d\n ", temp_word2, idngram_hash(temp_word2,M));*/

    add_to_idngram_hashtable(vocabulary,idngram_hash(temp_word2,M),temp_word2,vocab_size);
    if(vocab_size == M){
      quit(-1, "Number of entries reached the size of the hash.  Run the program again with a larger has size -hash \n");
    }
  }

  if (vocab_size > MAX_VOCAB_SIZE)    
    fprintf(stderr,"text2idngram : vocab_size %d\n is larger than %d\n",vocab_size,MAX_VOCAB_SIZE);

  return 0;
}

static int ng; /* Static variable to be used in qsort */

int compare_ngrams(const void *ngram1,
		   const void *ngram2
		   ) {

  int i;

  wordid_t *ngram_pointer1;
  wordid_t *ngram_pointer2;

  ngram_pointer1 = (wordid_t *) ngram1;
  ngram_pointer2 = (wordid_t *) ngram2;

  for (i=0;i<=ng-1;i++) {
    if (ngram_pointer1[i]<ngram_pointer2[i]) 
      return(-1);
    else {
      if (ngram_pointer1[i]>ngram_pointer2[i]) 
	return(1);
    }
  }

  return(0);

}

void add_to_buffer(wordid_t word_index,
		   int ypos,
		   int xpos, 
		   wordid_t *buffer) 
{
  buffer[(ng*ypos)+xpos] = word_index;
}

wordid_t buffer_contents(int ypos,
			  int xpos, 
			  wordid_t *buffer) 
{
  return (buffer[(ng*ypos)+xpos]);
}

int get_word( FILE *fp , char *word ) {

  /* read word from stream, checking for read errors and EOF */

  int nread;
  int rt_val;

  rt_val = 0;

  nread = fscanf( fp,MAX_WORD_LENGTH_FORMAT, word );

  if ( nread == 1 )
    rt_val = 1;
  else {
    if ( nread != EOF ) 
      quit(-1, "Error reading file" );
  }
  return(rt_val);
}

/*
  @return number_of_tempfiles
 */
int  read_txt2ngram_buffer(FILE* infp, 
			   struct idngram_hash_table *vocabulary, 
			   int32 verbosity,
			   wordid_t *buffer,
			   int buffer_size,
			   unsigned int n,
			   char* temp_file_root,
			   char* temp_file_ext,
			   FILE* temp_file
			   )
{
  /* Read text into buffer */
  char temp_word[MAX_WORD_LENGTH];
  int position_in_buffer;
  int number_of_tempfiles;
  unsigned int i,j;
  wordid_t *placeholder;
  wordid_t *temp_ngram;
  int temp_count;

#if 1
  int tmpval;
#endif

  temp_ngram  = (wordid_t *) rr_malloc(sizeof(wordid_t)*n);
  placeholder = (wordid_t *) rr_malloc(sizeof(wordid_t)*n);

  ng=n;

  position_in_buffer = 0;
  number_of_tempfiles = 0;

  //tk: looks like things may croak if the corpus has less than n words
  //not that such a corpus would be useful anyway
  for (i=0;i<=n-1;i++) {
    get_word(infp,temp_word);
    /*
        fprintf(stderr,"%s \n",temp_word);
	fprintf(stderr,"%d \n",index2(vocabulary,temp_word));
        fflush(stderr);
    */
    add_to_buffer(index2(vocabulary,temp_word),0,i,buffer);
  }

  while (!rr_feof(infp)) {
    /* Fill up the buffer */
    pc_message(verbosity,2,"Reading text into the n-gram buffer...\n");
    pc_message(verbosity,2,"20,000 n-grams processed for each \".\", 1,000,000 for each line.\n");

    while ((position_in_buffer<buffer_size) && (!rr_feof(infp))) {
      position_in_buffer++;
      show_idngram_nlines(position_in_buffer,verbosity);

      for (i=1;i<=n-1;i++) 
	add_to_buffer(buffer_contents(position_in_buffer-1,i,buffer),
		      position_in_buffer,i-1,buffer);
      
      if (get_word(infp,temp_word) == 1) {
      /*
	fprintf(stderr,"%s \n",temp_word);
	fprintf(stderr,"%d \n",index2(vocabulary,temp_word));
	fflush(stderr);
      */
	add_to_buffer(index2(vocabulary,temp_word),position_in_buffer,
		      n-1,buffer);
      }
    }

    for (i=0;i<=n-1;i++) 
      placeholder[i] = buffer_contents(position_in_buffer,i,buffer);

    /* Sort buffer */
    
    pc_message(verbosity,2,"\nSorting n-grams...\n");    
    
    qsort((void*) buffer,(size_t) position_in_buffer,n*sizeof(wordid_t),compare_ngrams);

    /* Output the buffer to temporary BINARY file */    
    number_of_tempfiles++;

    sprintf(temp_word,"%s/%hu%s",temp_file_root,
	    number_of_tempfiles,temp_file_ext);

    pc_message(verbosity,2,"Writing sorted n-grams to temporary file %s\n",
	       temp_word);

    temp_file = rr_oopen(temp_word);

    for (i=0;i<=n-1;i++) {
      temp_ngram[i] = buffer_contents(0,i,buffer);
#if MAX_VOCAB_SIZE < 65535
      /* This check is well-meaning but completely useless since
	 buffer_contents() can never return something greater than
	 MAX_VOCAB_SIZE (dhuggins@cs, 2006-03) */
      if (temp_ngram[i] > MAX_VOCAB_SIZE)
	quit(-1,"Invalid trigram in buffer.\nAborting");
#endif
    }
    temp_count = 1;

    for (i=1;i<=position_in_buffer;i++) {

      tmpval=compare_ngrams(temp_ngram,&buffer[i*n]);

      /*      for(k=0;k<=n-1;k++){
	fprintf(stderr, "tmpval: %d k %d, temp_ngram %d, &buffer[i*n] %d\n",tmpval, k, temp_ngram[k], (&buffer[i*n])[k]);
	}*/

      if (!compare_ngrams(temp_ngram,&buffer[i*n])) 
	temp_count++;
      else {
	/*	printf("Have been here?\n");*/
	for (j=0;j<=n-1;j++) {
	  rr_fwrite((char*) &temp_ngram[j],sizeof(wordid_t),1,
		    temp_file,"temporary n-gram ids");
	  temp_ngram[j] = buffer_contents(i,j,buffer);
	}
	rr_fwrite((char*)&temp_count,sizeof(int),1,temp_file,
		  "temporary n-gram counts");

	/*	for(j=0 ; j<=n-1;j++)
	  fprintf(stderr,"%d ",temp_ngram[j]);
	  fprintf(stderr,"%d\n",temp_count);*/

	temp_count = 1;
      }
    }
    
    rr_oclose(temp_file);

    for (i=0;i<=n-1;i++) 
      add_to_buffer(placeholder[i],0,i,buffer);

    position_in_buffer = 0;

  }

  return number_of_tempfiles;
}


void merge_tempfiles (int start_file, 
		      int end_file, 
		      char *temp_file_root,
		      char *temp_file_ext,
		      int max_files,
		      FILE *outfile,
		      int n,
		      int verbosity) {

  FILE *new_temp_file;
  char *new_temp_filename;
  
  FILE **temp_file;
  char **temp_filename;
  char **current_ngram;
  char smallest_ngram[1000];
  int *current_ngram_count;
  flag *finished;
  flag all_finished;
  int temp_count;
  char temp_word[500];
  int i,j;
  
  pc_message(verbosity,2,"Merging temp files %d through %d...\n", start_file,
 	  end_file);
   /*
    * If we try to do more than max_files, then merge into groups,
    * then merge groups recursively.
    */
    if (end_file-start_file+1 > max_files) {
       int new_start_file, new_end_file;
       int n_file_groups = 1 + (end_file-start_file)/max_files;
 
       fprintf(stderr, "%d files to do, in %d groups\n", end_file-start_file,
 	      n_file_groups);
 
       new_temp_filename = (char *) rr_malloc(300*sizeof(char));
 
       /*
        * These n_file_groups sets of files will be done in groups of
        * max_files batches each, as temp files numbered
        * end_file+1 ... end_file+n_file_groups,
        * and then these will be merged into the final result.
        */
 
       for (i = 0; i < n_file_groups; i++) {
 	  /* do files i*max_files through min((i+1)*max_files-1,end_file); */
 	  new_start_file = start_file + (i*max_files);
 	  new_end_file = start_file + ((i+1)*max_files) - 1;
 	  if (new_end_file > end_file) new_end_file = end_file;
 	  
 	  sprintf(new_temp_filename,
 		  "%s/%hu%s",
 		  temp_file_root,
 		  end_file+i+1,
 		  temp_file_ext);
 
 	  new_temp_file = rr_oopen(new_temp_filename);
 
 	  merge_tempfiles(new_start_file,
 			  new_end_file,
 			  temp_file_root,
			  temp_file_ext,
 			  max_files,
 			  new_temp_file,
 			  n,
			  verbosity);
 
 	  rr_iclose(new_temp_file);
 
       }
 
       merge_tempfiles(end_file+1,
		       end_file+n_file_groups,
		       temp_file_root,
		       temp_file_ext,
		       max_files,
		       outfile,
		       n,
		       verbosity);
 
       return;
    }
    
   /*
    * We know we are now doing <= max_files.
    */
 
   temp_file = (FILE **) rr_malloc((end_file+1)*sizeof(FILE *));
   temp_filename = (char **) rr_malloc((end_file+1)*sizeof(char *));
   for (i=start_file;i<=end_file;i++) {
     temp_filename[i] = (char *) rr_malloc(300*sizeof(char));
   }
   current_ngram = (char **) rr_malloc((end_file+1)*sizeof(char *));
   for (i=start_file;i<=end_file;i++) {
     current_ngram[i] = (char *) rr_malloc(1000*sizeof(char));
   }
   current_ngram_count = (int *) rr_malloc((end_file+1)*sizeof(int));
   finished = (flag *) rr_malloc(sizeof(flag)*(end_file+1));
  
   /* Open all the temp files for reading */
   for (i=start_file;i<=end_file;i++) {
     sprintf(temp_filename[i],"%s/%hu%s",
	     temp_file_root,i,temp_file_ext);
     temp_file[i] = rr_iopen(temp_filename[i]);
   }
 
   /* Now go through the files simultaneously, and write out the appropriate
      ngram counts to the output file. */
 
   for (i=start_file;i<=end_file;i++) {
     finished[i] = 0;
     if (!rr_feof(temp_file[i])) {
       for (j=0;j<=n-1;j++) {
 	if (fscanf(temp_file[i],"%s",temp_word) != 1) {
 	  if (!rr_feof(temp_file[i]))
 	    quit(-1,"Error reading temp file %s\n",temp_filename[i]);
 	}else {
 	  if (j==0)
 	    strcpy(current_ngram[i],temp_word);
  	  else {
 	    strcat(current_ngram[i]," ");
 	    strcat(current_ngram[i],temp_word);
  	  }
  	}
       }
       if (fscanf(temp_file[i],"%d",&current_ngram_count[i]) != 1) {
	 if (!rr_feof(temp_file[i]))
	   quit(-1,"Error reading temp file %s\n",temp_filename[i]);
       }
     }
   }
   
   all_finished = 0;
   
   while (!all_finished) {
  
     /* Find the smallest current ngram */
 
     strcpy(smallest_ngram,"");
 
     for (i=start_file;i<=end_file;i++) {
       if (!finished[i]) {
	 if (strcmp(smallest_ngram,current_ngram[i]) > 0 ||
	     (smallest_ngram[0] == '\0'))
	   strcpy(smallest_ngram,current_ngram[i]);
       }
     }
     
     /* For each of the files that are currently holding this ngram,
        add its count to the temporary count, and read in a new ngram
        from the files. */
  
     temp_count = 0;
 
     for (i=start_file;i<=end_file;i++) {
       if (!finished[i]) {
 	if (!strcmp(smallest_ngram,current_ngram[i])) {
 	  temp_count += current_ngram_count[i];
 	  if (!rr_feof(temp_file[i])) {
 	    for (j=0;j<=n-1;j++) {
 	      if (fscanf(temp_file[i],"%s",temp_word) != 1) {
 		if (!rr_feof(temp_file[i])) {
 		  quit(-1,"Error reading temp file %s\n",temp_filename[i]);
 		}
 	      }else {
 		if (j==0)
 		  strcpy(current_ngram[i],temp_word);
  		else {
 		  strcat(current_ngram[i]," ");
 		  strcat(current_ngram[i],temp_word);
  		}
  	      }
 	    }
 	    if (fscanf(temp_file[i],"%d",&current_ngram_count[i]) != 1) {
 	      if (!rr_feof(temp_file[i])) {
 		quit(-1,"Error reading temp file count %s\n",
 		     temp_filename[i]);
  	      }
  	    }
 	  }
 
 	  /*
 	   * PWP: Note that the fscanf may have changed the state of
 	   * temp_file[i], so we re-ask the question rather than just
 	   * doing an "else".
 	   */
 	  if (rr_feof(temp_file[i])) {
 	    finished[i] = 1;
 	    all_finished = 1;
 	    for (j=start_file;j<=end_file;j++) {
 	      if (!finished[j]) {
 		all_finished = 0;
  	      }
  	    }
  	  }
  	}
        }
      }
 
     /*
      * PWP: We cannot conditionalize this on (!all_finished) because
      * if we do we may have lost the very last count.  (Consider the
      * case when several files have ran out of data, but the last
      * couple have the last count in them.)
      */
     if (fprintf(outfile,"%s %d\n",smallest_ngram,temp_count) < 0) {
       quit(-1,"Write error encountered while attempting to merge temporary files.\nAborting, but keeping temporary files.\n");
     }
   }
 
   for (i=start_file;i<=end_file;i++) {
     rr_iclose(temp_file[i]);
     remove(temp_filename[i]);
   }
    
   free(temp_file);
   for (i=start_file;i<=end_file;i++) {
      free(temp_filename[i]);
    }
   free(temp_filename);  

   for (i=start_file;i<=end_file;i++) {
      free(current_ngram[i]);
   }
   free(current_ngram);

  free(current_ngram_count);
  free(finished);
}



/* Stuff that is common to both text2idngram and wngram2idngram */


static int n; /* Declare it globally, so doesn't need to be passed
		     as a parameter to compare_ngrams. which might
		     cause qsort to choke */
#include "../liblmest/ngram.h"
#include "../liblmest/stats.h"

int compare_ngrams3(const void *ngram1,
		   const void *ngram2
		   ) {

  int i;

  wordid_t *ngram_pointer1;
  wordid_t *ngram_pointer2;

  ngram_pointer1 = (wordid_t *) ngram1;
  ngram_pointer2 = (wordid_t *) ngram2;

  for (i=0;i<n;i++) {
    if (ngram_pointer1[i]<ngram_pointer2[i]) 
      return(1);
    else {
      if (ngram_pointer1[i]>ngram_pointer2[i]) 
	return(-1);
    }
  }

  return(0);

}



void merge_idngramfiles (int start_file, 
		      int end_file, 
		      char *temp_file_root,
		      char *temp_file_ext,
		      int max_files,
		      FILE *outfile,
		      flag write_ascii,
		      int fof_size,
		      int n_order) {
  FILE *new_temp_file;
  char temp_string[1000];
  char *new_temp_filename;
  
  FILE **temp_file;
  char **temp_filename;
  wordid_t **current_ngram;
  wordid_t *smallest_ngram;
  wordid_t *previous_ngram;

  int *current_ngram_count;
  flag *finished;
  flag all_finished;
  int temp_count;
  int i,j;
  flag first_ngram;
  fof_t **fof_array;
  ngram_sz_t *num_kgrams;
  int *ng_count;
  int pos_of_novelty;
  
  n = n_order;
  
  pos_of_novelty = n; /* Simply for warning-free compilation */
  num_kgrams = (ngram_sz_t *) rr_calloc(n-1,sizeof(ngram_sz_t));
  ng_count = (int *) rr_calloc(n-1,sizeof(int));
  first_ngram = 1;
  
  previous_ngram = (wordid_t *) rr_calloc(n,sizeof(wordid_t));
  temp_file = (FILE **) rr_malloc(sizeof(FILE *) * (end_file-start_file+1));
  temp_filename = (char **) rr_malloc(sizeof(char *) * 
				      (end_file-start_file+1));

  /* should change to 2d array*/
  current_ngram = (wordid_t **) rr_malloc(sizeof(wordid_t *) * 
						(end_file-start_file+1));
  for (i=0;i<=end_file-start_file;i++) 
    current_ngram[i] = (wordid_t *) rr_malloc(sizeof(wordid_t)*n);

  current_ngram_count = (int *) rr_malloc(sizeof(int)*(end_file-start_file+1));

  finished = (flag *) rr_malloc(sizeof(flag)*(end_file-start_file+1));
  smallest_ngram = (wordid_t *) rr_malloc(sizeof(wordid_t)*n);

  /* should change to 2d array*/
  fof_array = (fof_t **) rr_malloc(sizeof(fof_t *)*(n-1));
  for (i=0;i<=n-2;i++) 
    fof_array[i] = (fof_t *) rr_calloc(fof_size+1,sizeof(fof_t));

  if (end_file-start_file+1 > max_files) {
    sprintf(temp_string,"%s/%hu%s",temp_file_root,
	    end_file+1,temp_file_ext);
    new_temp_filename = salloc(temp_string);
    new_temp_file = rr_oopen(new_temp_filename);
    merge_tempfiles(start_file,start_file+max_files-1,
		    temp_file_root,temp_file_ext,max_files,
		    new_temp_file,write_ascii,0);
    merge_tempfiles(start_file+max_files,end_file+1,
		    temp_file_root,temp_file_ext,max_files,
		    outfile,write_ascii,0);
  }else {

    /* Open all the temp files for reading */
    for (i=0;i<=end_file-start_file;i++) {
      sprintf(temp_string,"%s/%hu%s",temp_file_root,
	      i+start_file,temp_file_ext);
      temp_filename[i] = salloc(temp_string);
      temp_file[i] = rr_iopen(temp_filename[i]);
    }
    
    /* Now go through the files simultaneously, and write out the appropriate
       ngram counts to the output file. */

    for (i=end_file-start_file;i>=0;i--) {
      finished[i] = 0;
      if (!rr_feof(temp_file[i])) {
	for (j=0;j<=n-1;j++) {
	  rr_fread((char*) &current_ngram[i][j], sizeof(wordid_t),1,
		   temp_file[i],"temporary n-gram ids",0);
	}    
	rr_fread((char*) &current_ngram_count[i], sizeof(int),1,
		 temp_file[i],"temporary n-gram counts",0);
      }
    }
    
    all_finished = 0;

    while (!all_finished) {

      /* Find the smallest current ngram */
      for (i=0;i<=n-1;i++) 
	smallest_ngram[i] = MAX_WORDID;

      for (i=0;i<=end_file-start_file;i++) {
	if (!finished[i]) {
	  if (compare_ngrams3(smallest_ngram,current_ngram[i]) < 0) {
	    for (j=0;j<n;j++)
	      smallest_ngram[j] = current_ngram[i][j];
	  }
	}
      }

#if MAX_VOCAB_SIZE < 65535
      /* This check is well-meaning but completely useless since
	 smallest_ngram[i] by definition cannot contain any value
	 greater than MAX_VOCAB_SIZE (dhuggins@cs, 2006-03) */
      for (i=0;i<=n-1;i++) {
	if (smallest_ngram[i] > MAX_VOCAB_SIZE) {
	  quit(-1,"Error : Temporary files corrupted, invalid n-gram found.\n");
	}
      }
#endif
	  
      /* For each of the files that are currently holding this ngram,
	 add its count to the temporary count, and read in a new ngram
	 from the files. */

      temp_count = 0;

      for (i=0;i<=end_file-start_file;i++) {
	if (!finished[i]) {
	  if (compare_ngrams3(smallest_ngram,current_ngram[i]) == 0) {
	    temp_count = temp_count + current_ngram_count[i];
	    if (!rr_feof(temp_file[i])) {
	      for (j=0;j<=n-1;j++) {
		rr_fread((char*) &current_ngram[i][j],sizeof(wordid_t),1,
			 temp_file[i],"temporary n-gram ids",0);
	      }
	      rr_fread((char*)&current_ngram_count[i],sizeof(int),1,
		       temp_file[i],"temporary n-gram count",0);
	    }else {
	      finished[i] = 1;
	      all_finished = 1;
	      for (j=0;j<=end_file-start_file;j++) {
		if (!finished[j]) 
		  all_finished = 0;
	      }
	    }
	  }
	}
      }
      
      if (write_ascii) {
	for (i=0;i<=n-1;i++) {

	  if (fprintf(outfile,"%d ",smallest_ngram[i]) < 0) 
	    {
	      quit(-1,"Write error encountered while attempting to merge temporary files.\nAborting, but keeping temporary files.\n");
	    }
	}
	if (fprintf(outfile,"%d\n",temp_count) < 0)  
	  quit(-1,"Write error encountered while attempting to merge temporary files.\nAborting, but keeping temporary files.\n");

      }else {
	for (i=0;i<=n-1;i++) {
	  rr_fwrite((char*)&smallest_ngram[i],sizeof(wordid_t),1,
		    outfile,"n-gram ids");
	}
	rr_fwrite((char*)&temp_count,sizeof(count_t),1,outfile,"n-gram counts");		   
      }

      if (fof_size > 0 && n>1) { /* Add stuff to fof arrays */
	
	/* Code from idngram2stats */	
	pos_of_novelty = n;
	for (i=0;i<=n-1;i++) {
	  if (smallest_ngram[i] > previous_ngram[i]) {
	    pos_of_novelty = i;
	    i=n;
	  }
	}
	  
	/* Add new N-gram */
	  
	num_kgrams[n-2]++;
	if (temp_count <= fof_size)
	  fof_array[n-2][temp_count]++;

	if (!first_ngram) {
	  for (i=n-2;i>=MAX(1,pos_of_novelty);i--) {
	    num_kgrams[i-1]++;
	    if (ng_count[i-1] <= fof_size) {
	      fof_array[i-1][ng_count[i-1]]++;
	    }
	    ng_count[i-1] = temp_count;
	  }
	}else {
	  for (i=n-2;i>=MAX(1,pos_of_novelty);i--) {
	    ng_count[i-1] = temp_count;
	  }
	  first_ngram = 0;
	}
	  
	for (i=0;i<=pos_of_novelty-2;i++)
	  ng_count[i] += temp_count;

	for (i=0;i<=n-1;i++)
	  previous_ngram[i]=smallest_ngram[i];

      }
    }
    
    for (i=0;i<=end_file-start_file;i++) {
      fclose(temp_file[i]);
      remove(temp_filename[i]); 
    }
    
  }    

  if (fof_size > 0 && n>1) { /* Display fof arrays */

    /* Process last ngram */
    for (i=n-2;i>=MAX(1,pos_of_novelty);i--) {
      num_kgrams[i-1]++;
      if (ng_count[i-1] <= fof_size)
	fof_array[i-1][ng_count[i-1]]++;

      ng_count[i-1] = temp_count;
    }
    
    for (i=0;i<=pos_of_novelty-2;i++)
      ng_count[i] += temp_count;

    display_fof_array(num_kgrams,fof_array,fof_size,stderr, n);

  }

}



