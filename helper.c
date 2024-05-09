#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <math.h>

static Py_ssize_t min(Py_ssize_t m, Py_ssize_t n) {
	return (m<n)?m:n;
}

static Py_ssize_t max(Py_ssize_t m, Py_ssize_t n) {
	return (m>n)?m:n;
}

static Py_ssize_t compute_alignment_score(const char* seq1, const char* seq2,
Py_ssize_t len1, Py_ssize_t len2, Py_ssize_t LEAD_GAPS) {
    if (len2 > len1) {
        /* swap the sequences so that seq1 is not shorter */
        const char* temp_seq;
        Py_ssize_t temp_len;
        temp_seq = seq1;
        seq1 = seq2;
        seq2 = temp_seq;
        temp_len = len1;
        len1 = len2;
        len2 = temp_len;
    }
    Py_ssize_t row[len2+1];
	Py_ssize_t pos1, pos2;
    for(pos2=0; pos2<=len2; pos2++) {
        /* initialize the row; allow LEAD_GAPS free leading gaps */
        row[pos2] = min(LEAD_GAPS-pos2, 0);
    }
    for(pos1=1; pos1<=len1; pos1++) {
        Py_ssize_t pre_score_diag = row[0];
        row[0] = min(LEAD_GAPS-pos1, 0);
        for(pos2=1; pos2<=len2; pos2++) {
            Py_ssize_t score_up = row[pos2-1] - (pos1==len1?0:1);
            Py_ssize_t score_diag = pre_score_diag + (seq1[pos1-1]==seq2[pos2-1]?1:-1);
            Py_ssize_t score_left = row[pos2] - (pos2==len2?0:1);
            pre_score_diag = row[pos2];
            row[pos2] = max(score_diag, max(score_up, score_left));
        }
    }
	/* check for sufficient identity */
    Py_ssize_t score = row[len2];
    return score;
}


//static PyObject *
//align_tails_v1(PyObject *self, PyObject *args)
//{
	///* initial input parsing and allocation */
	//Py_ssize_t num_args = PyTuple_Size(args);
	//if(num_args<2) {
        //if(!PyErr_Occurred()) 
            //PyErr_SetString(PyExc_TypeError,
            //"At least two sequences exptected");
        //return NULL;
    //}
	//const char** sequences = malloc(num_args * sizeof(char*));
	//if (!sequences) {
		//return PyErr_NoMemory();
	//}
	//Py_ssize_t* L = malloc(num_args * sizeof(Py_ssize_t));
	//if (!L) {
		//free(sequences);
		//return PyErr_NoMemory();
	//}
	///* pointers to the sequences */
	//Py_ssize_t n;
	//for(n=0; n<num_args; n++) {
		//PyObject* nth_arg = PyTuple_GetItem(args, n);
		//sequences[n] = PyBytes_AsString(nth_arg);
		//L[n] = strlen(sequences[n]);
	//}
	///* alignments */
	//Py_ssize_t numseq1=1, numseq2=1;
	//numseq1 += check_identity_c(sequences[0], sequences[1], L[0], L[1]);
	//for(n=2; n<num_args; n++) {
		//numseq1 += check_identity_c(sequences[0], sequences[n], L[0], L[n]);
		//numseq2 += check_identity_c(sequences[1], sequences[n], L[1], L[n]);
	//}
	//const char* seq = (numseq2>numseq1)?sequences[1]:sequences[0];
	//Py_ssize_t num_reads = max(numseq1, numseq2);
	///* final deallocation and return */
    //free(sequences);
    //free(L);
    //return Py_BuildValue("ny", num_reads, seq);
//}


static PyObject *
check_identity_c(PyObject *self, PyObject *args)
{
    char *seq1, *seq2;
    Py_ssize_t len1, len2, len_cutoff, lead_gaps;
    float pident;
    if (!PyArg_ParseTuple(args, "y#y#fnn", &seq1, &len1, &seq2, &len2,
    &pident, &len_cutoff, &lead_gaps)) {
        return NULL;
    }
    len1 = min(len1, len_cutoff);
    len2 = min(len2, len_cutoff);
    Py_ssize_t score = compute_alignment_score(seq1, seq2, len1, len2,
        lead_gaps);
    Py_ssize_t shorter_len = min(len1, len2);
    Py_ssize_t score_drop = shorter_len - score;
    double allowed_mism = (1 - pident) * shorter_len;
    if (score_drop <= 2 * allowed_mism) {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

static PyObject *
polya_score_c(PyObject *self, PyObject *args)
{
    int trim, s_match, s_mism, max_score=0, score=0;
    char *seq;
    Py_ssize_t len, pos;
    if (!PyArg_ParseTuple(args, "y#iii", &seq, &len, &trim, &s_match, &s_mism)) {
        return NULL;
    }
    for(pos=0; pos < len; pos++) {
        if (seq[pos] == 'T') {
            score += s_match;
        } else {
            score -= s_mism;
            if (pos < trim && score < 0) {
                score = 0;
            }
        }
        if (score > max_score) {
            max_score = score;
        }
    }
    return Py_BuildValue("i", max_score);
}

static PyObject *
polya_score_qual_c(PyObject *self, PyObject *args)
{
    short int *s_match, *s_mism;
    char *seq, *qual, *ch_s_match, *ch_s_mism;
    int trim, Q_OFFSET;
    Py_ssize_t len, qlen, len_ch_match, len_match, len_ch_mism, len_mism, max_q,
        pos, score=0, max_score=0;
    if (!PyArg_ParseTuple(args, "y#y#y#y#ii", &seq, &len, &qual, &qlen,
    &ch_s_match, &len_ch_match, &ch_s_mism, &len_ch_mism, &trim, &Q_OFFSET)) {
        return NULL;
    }
    s_match = (short int*) ch_s_match;
    s_mism = (short int*) ch_s_mism;
    len_match = len_ch_match * sizeof(char) / sizeof(short int);
    len_mism = len_ch_mism * sizeof(char) / sizeof(short int);
    max_q = min(len_match, len_mism) - 1;
    /* iterate through seq and qual */
    for(pos=0; pos < len; pos++) {
        Py_ssize_t q = min(max_q, qual[pos] - Q_OFFSET);
        if (seq[pos] == 'T') {
            score += s_match[q];
        } else {
            score -= s_mism[q];
            if (pos < trim && score < 0) {
                score = 0;
            }
        }
        if (score > max_score) {
            max_score = score;
        }
    }
    return Py_BuildValue("n", max_score);
}

static PyObject *
polya_score_full_length_c(PyObject *self, PyObject *args)
{
    int s_match, s_mism, score=0;
    char *seq;
    Py_ssize_t len, pos;
    if (!PyArg_ParseTuple(args, "y#ii", &seq, &len, &s_match, &s_mism)) {
        return NULL;
    }
    for(pos=0; pos < len; pos++) {
        if (seq[pos] == 'T') {
            score += s_match;
        } else {
            score -= s_mism;
        }
    }
    return Py_BuildValue("i", score);
}

static PyObject *
polya_score_full_length_qual_c(PyObject *self, PyObject *args)
{
    short int *s_match, *s_mism;
    char *seq, *qual, *ch_s_match, *ch_s_mism;
    int Q_OFFSET;
    Py_ssize_t len, qlen, len_ch_match, len_match, len_ch_mism, len_mism, max_q,
        pos, score=0;
    if (!PyArg_ParseTuple(args, "y#y#y#y#i", &seq, &len, &qual, &qlen,
    &ch_s_match, &len_ch_match, &ch_s_mism, &len_ch_mism, &Q_OFFSET)) {
        return NULL;
    }
    s_match = (short int*) ch_s_match;
    s_mism = (short int*) ch_s_mism;
    len_match = len_ch_match * sizeof(char) / sizeof(short int);
    len_mism = len_ch_mism * sizeof(char) / sizeof(short int);
    max_q = min(len_match, len_mism) - 1;
    /* iterate through seq and qual */
    for(pos=0; pos < len; pos++) {
        Py_ssize_t q = min(max_q, qual[pos] - Q_OFFSET);
        if (seq[pos] == 'T') {
            score += s_match[q];
        } else {
            score -= s_mism[q];
        }
    }
    return Py_BuildValue("n", score);
}


static PyObject *
compute_clipping_from_qual_c(PyObject *self, PyObject *args)
{
    unsigned char *qual, min_qual;
    Py_ssize_t len, len_r, len_l;
    Py_ssize_t pos, score, max_score, max_score_pos, clip_L, clip_R;
    if (!PyArg_ParseTuple(args, "y#Bnn", &qual, &len, &min_qual, &len_l,
    &len_r)) {
        return NULL;
    }
    /* right end */
    score = 0;
    max_score = 0;
    max_score_pos = len - len_r - 1;
    for(pos = len -len_r; pos < len; pos++) {
        if (qual[pos] >= min_qual) {
            score += 1;
            if (score > max_score) {
                max_score = score;
                max_score_pos = pos;
            }
        } else {
            score -= 2;
        }
    }
    clip_R = max_score_pos - len + len_r + 1;
    /* left end */
    score = 0;
    max_score = 0;
    max_score_pos = len_l;
    for(pos = len_l - 1; pos >= 0; pos--) {
        if (qual[pos] >= min_qual) {
            score += 1;
            if (score > max_score) {
                max_score = score;
                max_score_pos = pos;
            }
        } else {
            score -= 2;
        }
    }
    clip_L =  len_l - max_score_pos;
    /* return values */
    return Py_BuildValue("nn", clip_L, clip_R);
}


static PyObject *
test_c(PyObject *self, PyObject *args)
{
    /* TODO: what is the sizeof(elements) - C vs array.array typecodes */
    short int* score;
    char *seq;
    Py_ssize_t len_ch, len, pos;
    if (!PyArg_ParseTuple(args, "y#", &seq, &len_ch)) {
        return NULL;
    }
    score = (short int*) seq;
    len = len_ch * sizeof(char) / sizeof(short int);
    for(pos=0; pos < len; pos++) {
        PySys_FormatStderr("%zi %d \n", pos, (int) score[pos]);
    }
    return Py_BuildValue("iii", *score, score[1], score[2]);
}


static PyMethodDef AlignMethods[] = {
    {"check_identity",  check_identity_c, METH_VARARGS,
      "Check identity between sequences"},
    {"polya_score",  polya_score_c, METH_VARARGS, "Alignment score to poly(A)"},
    {"polya_score_qual",  polya_score_qual_c, METH_VARARGS,
      "Alignment score to poly(A) taking base qualities into account"},
    {"polya_score_full_length",  polya_score_full_length_c, METH_VARARGS,
      "Alignment score to poly(A), full length"},
    {"polya_score_full_length_qual",  polya_score_full_length_qual_c, METH_VARARGS,
      "Alignment score to poly(A), full length, accounting for base qualities"},
    {"compute_clipping_from_qual",  compute_clipping_from_qual_c, METH_VARARGS,
        "compute clipping coordinates excluding bases clipped for quality"},
    {"test",  test_c, METH_VARARGS,
        "test"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef alignmodule = {
    PyModuleDef_HEAD_INIT,
    "_TRHelper",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    AlignMethods
};

PyMODINIT_FUNC
PyInit__TRHelper(void)
{
    return PyModule_Create(&alignmodule);
}
