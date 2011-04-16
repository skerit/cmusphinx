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

#define DEFAULT_MAX_FILES 20

#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "../liblmest/toolkit.h"
#include "../libs/pc_general.h"
#include "../libs/general.h"
#include "../win32/compat.h"

int cmp_strings(const void *string1,const void *string2) {
  
  char *s1;
  char *s2;
  
  s1 = *((char **) string1);
  s2 = *((char **) string2);

  return (strcmp(s1,s2));

}

void help_message()
{
    fprintf(stderr,"text2wngram - Convert a text stream to a word n-gram stream.\n");
    fprintf(stderr,"Usage : text2wngram [ -n 3 ]\n");
    fprintf(stderr,"                    [ -chars %d ]\n",STD_MEM*7000000/11);
    fprintf(stderr,"                    [ -words %d ]\n",STD_MEM*1000000/11);
    fprintf(stderr,"                    [ -gzip | -compress ]\n");
    fprintf(stderr,"                    [ -verbosity 2 ]\n");
    fprintf(stderr,"                    < .text > .wngram\n");
}

int main (int argc, char **argv) {

  int n;
  int verbosity;
  int max_files;
  int max_words;
  int max_chars;

  int current_word;
  int current_char;
  int start_char;		/* start boundary (possibly > than 0) */

  int no_of_spaces;
  int pos_in_string;

  int i;
  char *current_string;
  char current_temp_filename[500];
  int current_file_number;
  FILE *temp_file;

  flag text_buffer_full;

  char *text_buffer;
  char **pointers;

  char current_ngram[500];
  int current_count;

  int counter;

  char temp_directory[1000];
  char *temp_file_ext;

  flag words_set;
  flag chars_set;

  /* Process command line */

  verbosity = pc_intarg(&argc, argv,"-verbosity",DEFAULT_VERBOSITY);
  pc_message(verbosity,2,"text2wngram\n");

  report_version(&argc,argv);

  if (pc_flagarg( &argc, argv,"-help")) {
    help_message();
    exit(1);
  }

  n = pc_intarg(&argc, argv,"-n",DEFAULT_N);

  /*  max_words = pc_intarg(&argc, argv,"-words",STD_MEM*1000000/11);
  max_chars = pc_intarg(&argc, argv,"-chars",STD_MEM*7000000/11); */

  max_words = pc_intarg(&argc, argv,"-words",-1);
  max_chars = pc_intarg(&argc, argv,"-chars",-1);

  if (max_words == -1) {
    words_set = 0;
    max_words = STD_MEM*1000000/11;
  }else
    words_set = 1;

  if (max_chars == -1) {
    chars_set = 0;
    max_chars = STD_MEM*7000000/11; 
  }else
    chars_set = 1;
  
  max_files = pc_intarg(&argc, argv,"-files",DEFAULT_MAX_FILES);

  if (pc_flagarg(&argc,argv,"-compress"))
    temp_file_ext = salloc(".Z");
  else {
    if (pc_flagarg(&argc,argv,"-gzip"))
      temp_file_ext = salloc(".gz");
    else
      temp_file_ext = salloc("");
  }

  strcpy(temp_directory, "cmuclmtk-XXXXXX");
  if (mkdtemp(temp_directory) == NULL) {
     quit(-1, "Failed to create temporary folder: %s\n", strerror(errno));
  }

  pc_report_unk_args(&argc,argv,verbosity);
 
  if (words_set && !chars_set)
    max_chars = max_words * 7;

  if (!words_set && chars_set)
    max_words = max_chars / 7;

  /* If the last charactor in the directory name isn't a / then add one. */
  
  pc_message(verbosity,2,"n = %d\n",n);
  pc_message(verbosity,2,"Number of words in buffer = %d\n",max_words);
  pc_message(verbosity,2,"Number of chars in buffer = %d\n",max_chars);
  pc_message(verbosity,2,"Max number of files open at once = %d\n",max_files);
  pc_message(verbosity,2,"Temporary directory = %s\n",temp_directory);

  /* Allocate memory for the buffers */

  text_buffer = (char *) rr_malloc(sizeof(char)*max_chars);
  pc_message(verbosity,2,"Allocated %d bytes to text buffer.\n",
	     sizeof(char)*max_chars);

  pointers = (char **) rr_malloc(sizeof(char *)*max_words);
  pc_message(verbosity,2,"Allocated %d bytes to pointer array.\n",
	     sizeof(char *)*max_words);

  current_file_number = 0;

  current_word = 1;
  start_char = 0;
  current_char = 0;
  counter = 0;
  pointers[0] = text_buffer;
      
  while (!feof(stdin)) {

    current_file_number++;

    /* Read text into buffer */
    
    pc_message(verbosity,2,"Reading text into buffer...\n");

    pc_message(verbosity,2,"Reading text into the n-gram buffer...\n");
    pc_message(verbosity,2,"20,000 words processed for each \".\", 1,000,000 for each line.\n");
    
    pointers[0] = text_buffer;
    
    while ((!rr_feof(stdin)) && 
	   (current_word < max_words) && 
	   (current_char < max_chars)) {

      text_buffer[current_char] = getchar();
      if (text_buffer[current_char] == '\n' || 
	  text_buffer[current_char] == '\t' ) {
	text_buffer[current_char] = ' ';
      }
      if (text_buffer[current_char] == ' ') {
	if (current_char > start_char) {
	  if (text_buffer[current_char-1] == ' ') {
	    current_word--;
	    current_char--;
	  }
	  pointers[current_word] = &(text_buffer[current_char+1]);
	  current_word++; 
	  counter++;
	  if (counter % 20000 == 0) {
	    if (counter % 1000000 == 0)
	      pc_message(verbosity,2,"\n");
	    else
	      pc_message(verbosity,2,".");
	  }
	}
      }
      
      if (text_buffer[current_char] != ' ' || current_char > start_char) 
	current_char++;
    }

    text_buffer[current_char]='\0';


    if (current_word == max_words || rr_feof(stdin)) {
      for (i=current_char+1;i<=max_chars-1;i++)
	text_buffer[i] = ' ';

      text_buffer_full = 0;
    }else
      text_buffer_full = 1;
    
    /* Sort buffer */

    pc_message(verbosity,2,"\nSorting pointer array...\n"); 

    qsort((void *) pointers,(size_t) current_word-n,sizeof(char *),cmp_strings);
   
    /* Write out temporary file */

    sprintf(current_temp_filename,"%s/%hu%s",temp_directory, current_file_number, temp_file_ext);

    pc_message(verbosity,2,"Writing out temporary file %s...\n",current_temp_filename);
        
    temp_file = rr_oopen(current_temp_filename);
    text_buffer[current_char] = ' ';
    
    current_count = 0;
    strcpy(current_ngram,"");
    
    for (i = 0; i <= current_word-n; i++) {
      current_string = pointers[i];
      
      /* Find the nth space */

      no_of_spaces = 0;
      pos_in_string = 0;
      while (no_of_spaces < n) {	
	if (current_string[pos_in_string] == ' ')
	  no_of_spaces++;

	pos_in_string++;
      }
      
      if (!strncmp(current_string,current_ngram,pos_in_string))
	current_count++;
      else {
	if (strcmp(current_ngram,""))
	  if (fprintf(temp_file,"%s %d\n",current_ngram,current_count) < 0) 
	    quit(-1,"Error writing to temporary file %s\n",current_temp_filename);

	current_count = 1;
	strncpy(current_ngram,current_string,pos_in_string);
	current_ngram[pos_in_string] = '\0';
      }
    }
    
    rr_oclose(temp_file);

    /* Move the last n-1 words to the beginning of the buffer, and set
       correct current_word and current_char things */

    strcpy(text_buffer,pointers[current_word-n]);
    pointers[0]=text_buffer;
   
    /* Find the (n-1)th space */

    no_of_spaces=0;
    pos_in_string=0;

    if (!text_buffer_full){ 
      while (no_of_spaces<(n-1)) {
	if (pointers[0][pos_in_string]==' ') {
	  no_of_spaces++;
	  pointers[no_of_spaces] = &pointers[0][pos_in_string+1];
	}
	pos_in_string++;
      }
    }else {
      while (no_of_spaces<n) {
	if (pointers[0][pos_in_string]==' ') {
	  no_of_spaces++;
	  pointers[no_of_spaces] = &pointers[0][pos_in_string+1];
	}
	pos_in_string++;
      }
      pos_in_string--;
    }

    current_char = pos_in_string;
    current_word = n;
    /* mark boundary beyond which counting pass cannot backup */
    start_char = current_char;

  }
  /* Merge temporary files */

  pc_message(verbosity,2,"Merging temporary files...\n");

  merge_tempfiles(1,
		  current_file_number,
		  temp_directory,
		  temp_file_ext,
		  max_files,
		  stdout,
		  n,
		  verbosity); 

  rmdir(temp_directory);
  pc_message(verbosity,0,"text2wngram : Done.\n");
  
  return 0;
}

