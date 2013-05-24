#include <Python.h>

static char module_docstring[] = "A docstring.";
static char align_docstring[] = "A docstring.";

static PyObject* nw_align(PyObject* self, PyObject* args) {
  // Read in arguments.
  char *seq1;
  char *seq2;
  int match;
  int mismatch;
  int gap;
  if (!PyArg_ParseTuple(args, "ssiii", &seq1, &seq2, &match, &mismatch, &gap))
    return NULL;

  int i;  // Used to iterate through sequence 1.
  int j;  // Used to iterate through sequence 2.
  int m = strlen(seq1) + 1;
  int n = strlen(seq2) + 1;

  // Declare arrays to be filled out.
  int maxAlignLen = m + n; // 1 too long, but makes life easier.
  char* align1 = (char*) malloc((maxAlignLen * sizeof(char)));
  char* align2 = (char*) malloc((maxAlignLen * sizeof(char)));
  int scores[n][m];
  int paths[n][m];

  // Set the first row and column of scores and paths.
  for (i=0; i<m; i++) {
    scores[0][i] = i * gap;
    paths[0][i] = 2; }
  for (j=0; j<n; j++) {
    scores[j][0] = j * gap;
    paths[j][0] = 3; }
  paths[0][0] = 0;

  /* Fill out scores and paths matrices. 
     In the paths matrix, 1 means diagonal, 2 is left, 3 is up, 0 is end. */
  int diag;
  int left;
  int up;
  for (j=1; j<n; j++) {
    for (i=1; i<m; i++) {
      if (seq1[i-1] == seq2[j-1])
	diag = scores[j-1][i-1] + match;
      else
	diag = scores[j-1][i-1] + mismatch;
      left = scores[j][i-1] + gap;
      up = scores[j-1][i] + gap;

      if (diag >= left && diag >= up) {
	scores[j][i] = diag;
	paths[j][i] = 1; }
      else if (left >= diag && left >= up) {
	scores[j][i] = left;
	paths[j][i] = 2; }
      else {
	scores[j][i] = up;
	paths[j][i] = 3; }
    }
  }

  /* DEBUG -- Print out matrix.
  for (j=1; j<n; j++) {
    for (i=1; i<m; i++) {
      printf("%i ", scores[j][i]); }
    printf("\n"); }
  */

  // Backtracking and filling out align arrays backwards.
  int alignPos = m + n - 1;
  align1[alignPos] = '\0';
  align2[alignPos] = '\0';
  alignPos--;
  i = m - 2;
  j = n - 2;
  int path = paths[j + 1][i + 1];
  while (path != 0) {
    if (path == 1) {
      align1[alignPos] = seq1[i];
      align2[alignPos] = seq2[j];
      i--;
      j--; }
    else if (path == 2) {
      align1[alignPos] = seq1[i];
      align2[alignPos] = '-';
      i--; }
    else {
      align1[alignPos] = '-';
      align2[alignPos] = seq2[j];
      j--; }
    alignPos--;
    path = paths[j + 1][i + 1];
  }

  // Copy the last x bytes from both align arrays to returnable objects.
  int alignLen = maxAlignLen - alignPos - 1;
  char* ret1 = (char*) malloc((alignLen * sizeof(char)));
  char* ret2 = (char*) malloc((alignLen * sizeof(char)));
  memcpy(ret1, align1 + (maxAlignLen-alignLen)*sizeof(char), alignLen*sizeof(char));
  memcpy(ret2, align2 + (maxAlignLen-alignLen)*sizeof(char), alignLen*sizeof(char));
  PyObject *ret = Py_BuildValue("(ss)", ret1, ret2);

  // Free the memory here.
  free(align1);
  free(align2);
  free(ret1);
  free(ret2);

  return ret;
}


static PyMethodDef module_methods[] = {
  {"align", nw_align, METH_VARARGS, align_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initnwmodule(void) {
  PyObject *m = Py_InitModule3("nwmodule", module_methods, module_docstring);
  if (m == NULL)
    return;
}
