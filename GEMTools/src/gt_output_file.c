/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_file.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include <sys/wait.h>
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif
#include "gt_output_file.h"
#include "gt_pipe_io.h"

/*
 * Setup
 */
GT_INLINE void gt_output_file_init_buffers(gt_output_file* const output_file) {
  /* Output Buffers */
  uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS;++i) {
    output_file->buffer[i]=NULL;
  }
  output_file->buffer_busy=0;
  output_file->buffer_write_pending=0;
  /* Block ID (for synchronization purposes) */
  output_file->mayor_block_id=0;
  output_file->minor_block_id=0;
  /* Mutexes */
  gt_cond_fatal_error(pthread_cond_init(&output_file->out_buffer_cond,NULL),SYS_COND_VAR_INIT);
  gt_cond_fatal_error(pthread_cond_init(&output_file->out_write_cond,NULL),SYS_COND_VAR_INIT);
  gt_cond_fatal_error(pthread_mutex_init(&output_file->out_file_mutex, NULL),SYS_MUTEX_INIT);
}

#ifdef HAVE_ZLIB
static void* gt_output_file_pipe_gzip(void *s)
{
	u_int8_t *buffer[GT_OUTPUT_COMPRESS_BUFFER_SIZE];
	gt_output_file* output_file=s;
	FILE *in;
  gt_cond_fatal_error(!(in=fdopen(output_file->pipe_fd[0],"r")),FILE_FDOPEN);
  while(!feof(in)) {
  	size_t n=fread(buffer,1,GT_OUTPUT_COMPRESS_BUFFER_SIZE,in);
  	if(!n) break;
  	size_t bytes_written=gzwrite(output_file->cfile,buffer,n);
    gt_cond_fatal_error(bytes_written!=n,OUTPUT_FILE_FAIL_WRITE);
  }
  int err=gzclose(output_file->cfile);
  gt_cond_error(err!=Z_OK,FILE_CLOSE,output_file->file_name);
  fclose(in);
	return 0;
}
#endif

#ifdef HAVE_BZLIB
static void* gt_output_file_pipe_bzip(void *s)
{
	u_int8_t *buffer[GT_OUTPUT_COMPRESS_BUFFER_SIZE];
	gt_output_file* output_file=s;
	FILE *in;
  gt_cond_fatal_error(!(in=fdopen(output_file->pipe_fd[0],"r")),FILE_FDOPEN);
  while(!feof(in)) {
  	size_t n=fread(buffer,1,GT_OUTPUT_COMPRESS_BUFFER_SIZE,in);
  	if(!n) break;
  	int err;
  	BZ2_bzWrite(&err,output_file->cfile,buffer,n);
    gt_cond_fatal_error(err!=BZ_OK,OUTPUT_FILE_FAIL_WRITE);
  }
  int err;
  BZ2_bzWriteClose(&err,output_file->cfile,0,NULL,NULL);
  gt_cond_error(err!=BZ_OK,FILE_CLOSE,output_file->file_name);
  fclose(in);
	return 0;
}
#endif

gt_output_file* gt_output_stream_new_compress(FILE* const file,const gt_output_file_type output_file_type, gt_output_file_compression compression_type) {

  GT_NULL_CHECK(file);
  gt_output_file* output_file = gt_alloc(gt_output_file);
  /* Output file */
  output_file->file_name=GT_STREAM_FILE_NAME;
  output_file->file=file;
  output_file->file_type=output_file_type;

#ifndef HAVE_ZLIB
  if(compression_type==GZIP) compression_type=NONE;
#endif
#ifndef HAVE_BZLIB
  if(compression_type==BZIP2) compression_type=NONE;
#endif
  if(compression_type!=NONE && isatty(fileno(file))) {
  	fprintf(stderr,"Will not output compressed data to a tty\n");
  	compression_type=NONE;
  }
  output_file->compression_type=compression_type;
  switch(compression_type) {
  case GZIP:
#ifdef HAVE_ZLIB
		gt_cond_fatal_error(pipe(output_file->pipe_fd)<0,SYS_PIPE);
 	  gt_cond_fatal_error(!(output_file->cfile=gzdopen(dup(fileno(file)),"wb")),FILE_GZIP_OPEN,output_file->file_name);
	  gt_cond_fatal_error(pthread_create(&output_file->pth,NULL,gt_output_file_pipe_gzip,output_file),SYS_THREAD);
	  gt_cond_fatal_error(!(output_file->file=fdopen(output_file->pipe_fd[1],"w")),FILE_FDOPEN);
#endif
	  break;
  case BZIP2:
#ifdef HAVE_BZLIB
		gt_cond_fatal_error(pipe(output_file->pipe_fd)<0,SYS_PIPE);
	  int er;
	  output_file->cfile=BZ2_bzWriteOpen(&er,output_file->file,1,0,0);
	  gt_cond_fatal_error(er!=BZ_OK,FILE_BZIP2_OPEN,output_file->file_name);
	  gt_cond_fatal_error(pthread_create(&output_file->pth,NULL,gt_output_file_pipe_bzip,output_file),SYS_THREAD);
	  gt_cond_fatal_error(!(output_file->file=fdopen(output_file->pipe_fd[1],"w")),FILE_FDOPEN);
#endif
  	break;
  default:
  	break;
  }
  /* Setup buffers */
  gt_output_file_init_buffers(output_file);
  return output_file;
}

gt_output_file *gt_output_file_pipe_open(char *const file_name,const gt_output_file_type output_file_type)
{
  char *trimmed_name=NULL;
  gt_cond_fatal_error(gt_pipe_io_check_command(file_name,&trimmed_name)!=GT_PIPE_IO_WRITE,PIPE_NOT_WRITE,file_name);
  GT_NULL_CHECK(trimmed_name);
  int fd=gt_pipe_io_child_open(GT_PIPE_IO_WRITE,trimmed_name);
  gt_cond_fatal_error(fd<0,SYS_PIPE); // Should never happen
  FILE *stream=fdopen(fd,"w");
  gt_cond_fatal_error(stream==NULL,FILE_OPEN,trimmed_name); // Should never happen
  free(trimmed_name);
  return gt_output_stream_new_compress(stream,output_file_type,NONE); // We ignore compression if we writing to a pipe
}

gt_output_file* gt_output_file_new_compress(char* const file_name,const gt_output_file_type output_file_type,gt_output_file_compression compression_type)
{
  GT_NULL_CHECK(file_name);
  // Check for pipe
  if(strchr(file_name,'|')) return gt_output_file_pipe_open(file_name,output_file_type);
  // Allocate handler
  gt_output_file* output_file = gt_alloc(gt_output_file);
  // Output file
  output_file->file_name=file_name;
#ifndef HAVE_ZLIB
  if(compression_type==GZIP) compression_type=NONE;
#endif
#ifndef HAVE_BZLIB
  if(compression_type==BZIP2) compression_type=NONE;
#endif
	output_file->compression_type=compression_type;
  switch(compression_type) {
  	case GZIP:
#ifdef HAVE_ZLIB
  		gt_cond_fatal_error(pipe(output_file->pipe_fd)<0,SYS_PIPE);
  	  gt_cond_fatal_error(!(output_file->cfile=gzopen(file_name,"wb")),FILE_GZIP_OPEN,file_name);
  	  gt_cond_fatal_error(pthread_create(&output_file->pth,NULL,gt_output_file_pipe_gzip,output_file),SYS_THREAD);
  	  gt_cond_fatal_error(!(output_file->file=fdopen(output_file->pipe_fd[1],"w")),FILE_FDOPEN);
#endif
  	  break;
  	case BZIP2:
#ifdef HAVE_BZLIB
  		gt_cond_fatal_error(pipe(output_file->pipe_fd)<0,SYS_PIPE);
  	  gt_cond_fatal_error(!(output_file->file=fopen(file_name,"w")),FILE_OPEN,file_name);
  	  int er;
  	  output_file->cfile=BZ2_bzWriteOpen(&er,output_file->file,1,0,0);
  	  gt_cond_fatal_error(er!=BZ_OK,FILE_BZIP2_OPEN,file_name);
  	  gt_cond_fatal_error(pthread_create(&output_file->pth,NULL,gt_output_file_pipe_bzip,output_file),SYS_THREAD);
  	  gt_cond_fatal_error(!(output_file->file=fdopen(output_file->pipe_fd[1],"w")),FILE_FDOPEN);
#endif
  	  break;
  	default:
  		gt_cond_fatal_error(!(output_file->file=fopen(file_name,"w")),FILE_OPEN,file_name);
  		break;
  }
  output_file->file_type=output_file_type;
  /* Setup buffers */
  gt_output_file_init_buffers(output_file);
  return output_file;
}

gt_status gt_output_file_close(gt_output_file* const output_file) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  gt_status error_code = 0;
  switch(output_file->compression_type) {
  case GZIP:
  case BZIP2:
  	error_code|=fclose(output_file->file);
	  gt_cond_error(error_code,FILE_CLOSE,output_file->file_name);
  	pthread_join(output_file->pth,NULL);
  	break;
  default:
    // Close file not stream
    if(strcmp(output_file->file_name, GT_STREAM_FILE_NAME)) {
    	error_code|=fclose(output_file->file);
  	  gt_cond_error(error_code,FILE_CLOSE,output_file->file_name);
    } else {
      // In case the stream is from a pipe where we forked a shell we should wait for the child to die
      signal(SIGCHLD,SIG_DFL);
      int i;
      while(waitpid(-1,&i,WNOHANG)>0);
    }
  	break;
  }
  // Delete allocated buffers
  uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&output_file->buffer[i]!=NULL;++i) {
    gt_output_buffer_delete(output_file->buffer[i]);
  }
  // Free mutex/CV
  gt_cond_error(error_code|=pthread_cond_destroy(&output_file->out_buffer_cond),SYS_COND_VAR_INIT);
  gt_cond_error(error_code|=pthread_cond_destroy(&output_file->out_write_cond),SYS_COND_VAR_INIT);
  gt_cond_error(error_code|=pthread_mutex_destroy(&output_file->out_file_mutex),SYS_MUTEX_DESTROY);
  // Free handler
  gt_free(output_file);
  return error_code;
}

/*
 * Output File Printers
 */
GT_INLINE gt_status gt_vofprintf(gt_output_file* const output_file,const char *template,va_list v_args) {
  GT_OUTPUT_FILE_CHECK(output_file);
  GT_NULL_CHECK(template);
  gt_status error_code;
  GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
  {
    error_code = vfprintf(output_file->file,template,v_args);
  }
  GT_END_MUTEX_SECTION(output_file->out_file_mutex);
  return error_code;
}
GT_INLINE gt_status gt_ofprintf(gt_output_file* const output_file,const char *template,...) {
  GT_OUTPUT_FILE_CHECK(output_file);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  const gt_status error_code = gt_vofprintf(output_file,template,v_args);
  va_end(v_args);
  return error_code;
}

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* __gt_buffered_output_file_request_buffer(gt_output_file* const output_file) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  // Conditional guard. Wait till there is any free buffer left
  while (output_file->buffer_busy==GT_MAX_OUTPUT_BUFFERS) {
    GT_CV_WAIT(output_file->out_buffer_cond,output_file->out_file_mutex);
  }
  // There is at least one free buffer. Get it!
  uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&output_file->buffer[i]!=NULL;++i) {
    if (gt_output_buffer_get_state(output_file->buffer[i])==GT_OUTPUT_BUFFER_FREE) break;
  }
  gt_cond_fatal_error(i>=GT_MAX_OUTPUT_BUFFERS,ALG_INCONSISNTENCY);
  if (output_file->buffer[i]==NULL) {
    output_file->buffer[i] = gt_output_buffer_new();
  }
  ++output_file->buffer_busy;
  gt_output_buffer_initiallize(output_file->buffer[i],GT_OUTPUT_BUFFER_BUSY);
  return output_file->buffer[i];
}
GT_INLINE gt_output_buffer* gt_output_file_request_buffer(gt_output_file* const output_file) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  gt_output_buffer* fresh_buffer;
  GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
  {
    fresh_buffer = __gt_buffered_output_file_request_buffer(output_file);
  }
  GT_END_MUTEX_SECTION(output_file->out_file_mutex);
  return fresh_buffer;
}

GT_INLINE void __gt_buffered_output_file_release_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  // Broadcast. Wake up sleepy.
  if (output_file->buffer_busy==GT_MAX_OUTPUT_BUFFERS) {
    GT_CV_BROADCAST(output_file->out_buffer_cond);
  }
  // Free buffer
  gt_output_buffer_set_state(output_buffer,GT_OUTPUT_BUFFER_FREE);
  --output_file->buffer_busy;
  gt_output_buffer_clear(output_buffer);
}
GT_INLINE void gt_output_file_release_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
  {
    __gt_buffered_output_file_release_buffer(output_file,output_buffer);
  }
  GT_END_MUTEX_SECTION(output_file->out_file_mutex);
}

GT_INLINE gt_output_buffer* gt_output_file_write_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  if (gt_output_buffer_get_used(output_buffer) > 0) {
    int64_t bytes_written;
    gt_vector* const vbuffer = gt_output_buffer_to_vchar(output_buffer);
    GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
    {
      bytes_written = fwrite(gt_vector_get_mem(vbuffer,char),1,
          gt_vector_get_used(vbuffer),output_file->file);
    }
    GT_END_MUTEX_SECTION(output_file->out_file_mutex);
    gt_cond_fatal_error(bytes_written!=gt_vector_get_used(vbuffer),OUTPUT_FILE_FAIL_WRITE);
  }
  gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_BUSY);
  return output_buffer;
}
GT_INLINE gt_output_buffer* gt_output_file_sorted_write_buffer_asynchronous(
    gt_output_file* const output_file,gt_output_buffer* output_buffer,const bool asynchronous) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  // Set the block buffer as write pending and set the victim
  bool victim;
  uint32_t mayor_block_id, minor_block_id;
  GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
  {
    victim = (output_file->mayor_block_id==gt_output_buffer_get_mayor_block_id(output_buffer) &&
              output_file->minor_block_id==gt_output_buffer_get_minor_block_id(output_buffer));
    while (!asynchronous && !victim) {
      GT_CV_WAIT(output_file->out_write_cond,output_file->out_file_mutex);
      victim = (output_file->mayor_block_id==gt_output_buffer_get_mayor_block_id(output_buffer) &&
                output_file->minor_block_id==gt_output_buffer_get_minor_block_id(output_buffer));
    }
    // Set the buffer as write pending
    ++output_file->buffer_write_pending;
    gt_output_buffer_set_state(output_buffer,GT_OUTPUT_BUFFER_WRITE_PENDING);
    if (!victim) { // Enqueue the buffer and continue (someone else will do the writing)
      output_buffer = __gt_buffered_output_file_request_buffer(output_file);
      GT_END_MUTEX_SECTION(output_file->out_file_mutex);
      return output_buffer;
    } else {
      mayor_block_id = output_file->mayor_block_id;
      minor_block_id = output_file->minor_block_id;
    }
  }
  GT_END_MUTEX_SECTION(output_file->out_file_mutex);
  // I'm the victim, I will output as much as I can
  do {
    // Write the current buffer
    if (gt_output_buffer_get_used(output_buffer) > 0) {
      gt_vector* const vbuffer = gt_output_buffer_to_vchar(output_buffer);
      const int64_t bytes_written =
          fwrite(gt_vector_get_mem(vbuffer,char),1,gt_vector_get_used(vbuffer),output_file->file);
      gt_cond_fatal_error(bytes_written!=gt_vector_get_used(vbuffer),OUTPUT_FILE_FAIL_WRITE);
    }
    // Update buffers' state
    GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex) {
      // Decrement write-pending blocks and update next block ID (mayorID,minorID)
      --output_file->buffer_write_pending;
      if (output_buffer->is_final_block) {
        ++mayor_block_id;
        minor_block_id = 0;
      } else {
        ++minor_block_id;
      }
      // Search for the next block buffer in order
      bool next_found = false;
      if (output_file->buffer_write_pending>0) {
        uint64_t i;
        for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&output_file->buffer[i]!=NULL;++i) {
          if (mayor_block_id==gt_output_buffer_get_mayor_block_id(output_file->buffer[i]) &&
              minor_block_id==gt_output_buffer_get_minor_block_id(output_file->buffer[i])) {
            if (gt_output_buffer_get_state(output_file->buffer[i])!=GT_OUTPUT_BUFFER_WRITE_PENDING) {
              break; // Cannot dump a busy buffer
            } else {
              // I'm still the victim, free the current buffer and output the new one
              __gt_buffered_output_file_release_buffer(output_file,output_buffer);
              next_found = true;
              output_buffer = output_file->buffer[i];
              break;
            }
          }
        }
      }
      if (!next_found) {
        // Fine, I'm done, let's get out of here ASAP
        output_file->mayor_block_id = mayor_block_id;
        output_file->minor_block_id = minor_block_id;
        gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_BUSY);
        GT_CV_BROADCAST(output_file->out_write_cond);
        GT_END_MUTEX_SECTION(output_file->out_file_mutex);
        return output_buffer;
      }
    } GT_END_MUTEX_SECTION(output_file->out_file_mutex);
  } while (true);
}
GT_INLINE gt_output_buffer* gt_output_file_dump_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer,const bool asynchronous) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  switch (output_file->file_type) {
    case SORTED_FILE:
      return gt_output_file_sorted_write_buffer_asynchronous(output_file,output_buffer,asynchronous);
      break;
    case UNSORTED_FILE:
      return gt_output_file_write_buffer(output_file,output_buffer);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      break;
  }
}
