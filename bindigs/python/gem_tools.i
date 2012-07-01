%module gemtools
%{
#include "gem_tools.h"
#include <sys/stat.h>
%}

%include "typemaps.i"

// convert python files
%typemap(python,in) FILE * {
  if (!PyFile_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a file!");
      return NULL;
  }
  $1 = PyFile_AsFile($input);
}


typedef int gt_status;
typedef unsigned long uint64_t;


%include "gt_commons.h"

// Input handlers
%include "gt_input_file.h"
%include "gt_buffered_map_input.h"

// Output handlers
%include "gt_buffered_output_file.h"
//#include "gt_output_map.h"

// GEM-Tools basic data structures: Template/Alignment/Maps/...
%include "gt_misms.h"
%include "gt_map.h"
%include "gt_alignment.h"
%include "gt_template.h"

%include "gt_template.h"

