/** @file qcd_utils.h
 *  @brief header file for qcd_utils.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 *  @todo: rename get* --> read*
 */

#ifndef H_QCD_UTILS
#define H_QCD_UTILS 1 

/* prototypes */

void *qcd_malloc(size_t size, int p);
void qcd_free(void *addr);

/**
 * @fn int qcd_isBigEndian(void)
 * @brief Determines whether a machine is little or big endian.
 * @details A 1 is returned if the machine is big endian, a
 *          0 if it is little endian.
 * @author Giannis Koutsou, minor adjustments Tomasz Korzec
 * @return 0 or 1
 */
int qcd_isBigEndian(void);

/**
 * @fn void qcd_swap8(double *Rd, int N)
 * @brief Changes the byte ordering in a vector of doubles.
 * @details The N double precision numbers pointed to by Rd are changed
 *          from little to big endian or vice versa.
 * @param[in,out] Rd pointer to an array of N double precision floating 
 *                point numbers
 * @param[in] N   Length of the array Rd
 */
void qcd_swap8(double *Rd, int N);

/**
 * @fn void qcd_swap4(float *Rd, int N)
 * @brief Changes the byte ordering in a vector of floats.
 * @details The N single precision numbers pointed to by Rd are changed
 *          from little to big endian or vice versa.
 * @param[in,out] Rd pointer to an array of N single precision floating 
 *                point numbers
 * @param[in] N   Length of the array Rd
 */
void qcd_swap4(float *Rd, int N);

/**
 * @fn void qcd_swapInt4(int *dat, int N)
 * @brief Changes the byte ordering in a vector of integers.
 * @details The N integers pointed to by dat are changed
 *          from little to big endian or vice versa.
 * @param[in,out] dat pointer to an array of N integers (4 bytes)
 * @param[in] N   Length of the array dat
 */
void qcd_swapInt4(int *dat, int N);

/**
 * @fn char* qcd_getParam(char token[], char* params, long len)
 * @brief Looks for a particular parameter and returns its value.
 * @details The string params of length len is scanned for the first occurrence
 *          of the 0-terminated string token. The string starting at the first
 *          character after token until the next occurrence of the character '<' 
 *          is returned (as a 0-terminated string).
 * @param[in] token a string that marks the beginning of the parameter field
 * @param[in] params the complete input file as a string
 * @param[in] len the length of the string in params
 * @return  the parameter value, i.e. the string between token and the next '<'
 * @author Giannis Koutsou, minor adjustments Tomasz Korzec
 */
char* qcd_getParam(char token[], char* params, long len);

/**
 * @fn int qcd_findOpt(int argc,char *argv[],char *opt)
 * @brief Returns the position of the argument that matches the string opt
 * @details On process 0, this function compares the string opt with the 
 *          arguments argv[1],..,argv[argc-1] and returns the position of the 
 *          first argument that matches the string. If there is no matching 
 *          argument, or if the function is called from another process, the 
 *          return value is 0.
 * @param[in] argc number of arguments, as passed to main
 * @param[in] argv array of strings as passed to main
 * @param[in] opt string of the argument
 * @return  the position at which opt occurs in argv
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
int qcd_findOpt(int argc,char *argv[],char *opt);

/**
 * @fn int qcd_isDirWritable(char* dir)
 * @brief Tests whether writable files can be created in a directory
 * @details It is tested, if writable files can be opened in the directory
 *          dir. If yes 1 is returned, if no, 0.
 * @param[in] char string of the directory to test
 * @return 0 if no writable files can be created in dir, 1 otherwise
 * @author Martin Luescher, adapted to qcdlib: Tomasz Korzec
 */
int qcd_isDirWritable(char* dir);

/**
 * @fn int qcd_nameSize(char *format,...)
 * @brief Tests whether writable files can be created in a directory
 * @details On process 0, this function returns the length of the string that
*     would be printed by calling sprintf(*,format,...). The format
*     string can be any combination of literal text and the conversion
*     specifiers %s, %d and %.nf (where n is a positive integer). When
*     called on other processes, the function does nothing and returns
*     the value of qcd_MAX_STRING_LENGTH. In case of an error -1 is returned.
 * @param[in] format a string with formating information
 * @return the expected length of an output string or -1 in case of errors
 * @author Martin Luescher, adapted to qcdlib: Tomasz Korzec
 */
int qcd_nameSize(char *format,...);

/**
 * @fn long qcd_findOpenQCDSection(const char * title)
 * @brief returns the offset to the first line after the string "[title]"
 * @details On process 0, this function scans stdin for a line starting with
*     the string "[title]" (after any number of blanks). -1 is returned and an
 *    error message printed if no such line is found or if there are
*     several of them. The function returns the offset of the line from
*     the beginning of the file and positions the file pointer to the
*     next line. On processes other than 0, the function does nothing 
*     and returns -1.
 * @param[in] title string with the section name
 * @return offset of the line
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
long qcd_findOpenQCDSection(const char *title);


/**
 * @fn long qcd_readOpenQCDLine(const char *tag,const char *format,...)
 * @brief Reads data from a line with a tag
 * @details On process 0, this function reads a line of text and data from stdin
 *     in a controlled manner. The tag can
 *     be the empty string "" and must otherwise be an alpha-numeric word
 *     that starts with a letter. If it is not empty, the program searches
 *     for the tag in the current section. An error occurs if the tag is not
 *     found. The program returns the offset of the line from the beginning
 *     of the file and positions the file pointer to the next line. On
 *     processes other than 0, the program does nothing and returns -1L. In 
 *    case of an error also -1L is returned.
 * @param[in] tag string containing the tag
 * @param[in] format the format string that should match the line from the tag
 *            on.
 * @param[in,out] ... pointers to variables into which the data is written
 * @return offset of the line with the matched tag. -1 in case of an error
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
long qcd_readOpenQCDLine(const char *tag,const char *format,...);


/**
 * @fn long qcd_readOpenQCDLineOpt(const char *tag,const char *defaultline,const char *format, ...)
 * @brief Reads optional data from a line with a tag
 * @details On process 0, this program reads a line of text and data from stdin
 *     in a controlled manner. The tag cannot
 *     be the empty string "", the program searches
 *     for the tag in the current section. If the tag is not
 *     found, default values are assigned to the variables. The program returns
 *     the offset of the line from the beginning
 *     of the file and positions the file pointer to the next line. On
 *     processes other than 0, the program does nothing and returns -1L.
 * @param[in] tag string containing the tag
 * @param[in] defaultline a string from which the data is parsed in case the tag
 *            is not found in stdin.
 * @param[in] format the format string that should match the line from the tag
 *            on.
 * @param[in,out] ... pointers to variables into which the data is written
 * @return offset of the line with the matched tag. -1 in case of an error
 * @author Martin Luescher, adjustments Tomasz Korzec
 */
long qcd_readOpenQCDLineOpt(const char *tag,const char *defaultline,const char *format, ...);


/**
 * @fn int qcd_countOpenQCDtokens(const char *tag)
 * @brief Counts the nomber of tokens in a line with a tag
 * @details On process 0, this function finds and reads a line from stdin, exactly
 *     as qcd_readOpenQCDLine(tag,..) does, and returns the number of tokens found on
 *     that line after the tag. Tokens are separated by white space (blanks,
 *     tabs or newline characters) and comments (text beginning with #) are
 *     ignored. On exit, the file pointer is positioned at the next line. If
 *     called on other processes, the program does nothing and returns 0. 
 *     In case of errors -1 is returned.
 * @param[in] tag string containing the tag
 * @return an integer with the number of tokens or -1 in case of an error
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
int qcd_countOpenQCDtokens(const char *tag);

/**
 * @fn int qcd_readOpenQCDIprms(const char *tag, int n, int *iprms)
 * @brief Reads integer data from a line with a tag
 * @details On process 0, this function finds and reads a line from stdin, 
 *    exactly as qcd_readOpenQCDLine(tag,..) does, reads n integer values from 
 *    that line after the tag and assigns them to the elements of the array 
 *    iprms. An error occurs if less than n values are found on the line. The
 *    values must be separated by white space (blanks, tabs or newline 
 *    characters). On exit, the file pointer is positioned at the next line. 
 *    When called on other processes, the program does nothing. In case of
 *    an error 1 is returned, in case of success 0.
 * @param[in] tag string containing the tag
 * @param[in] n number of integer data items to read
 * @param[in,out] iprms pointer to an array of integers into which the read data
 *                is written.
 * @return 0 on success, 1 on failure
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
int qcd_readOpenQCDIprms(const char *tag, int n, int *iprms);


/**
 * @fn int qcd_readOpenQCDDprms(const char *tag, int n, int *dprms)
 * @brief Reads double precision data from a line with a tag
 * @details On process 0, this function finds and reads a line from stdin, 
 *    exactly as qcd_readOpenQCDLine(tag,..) does, reads n double values from 
 *    that line after the tag and assigns them to the elements of the array 
 *    dprms. An error occurs if less than n values are found on the line. The
 *    values must be separated by white space (blanks, tabs or newline 
 *    characters). On exit, the file pointer is positioned at the next line. 
 *    When called on other processes, the program does nothing. In case of
 *    an error 1 is returned, in case of success 0.
 * @param[in] tag string containing the tag
 * @param[in] n number of double data items to read
 * @param[in,out] dprms pointer to an array of double into which the read data
 *                is written.
 * @return 0 on success, 1 on failure
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
int qcd_readOpenQCDDprms(const char *tag, int n, double *dprms);

/**
 * @fn int qcd_copyFile(char *in, char *out)
 * @brief copies a file
 * @param[in] in file name of intput file
 * @param[in] out file name of output file
 * @return 0 on success, 1 on error
 */
int qcd_copyFile(char *in, char *out);

/**
 * @fn void qcd_printMemInfo()
 * @brief prints out a short memory usage statistics to stdout
 */
void qcd_printMemInfo(void);

/**
 * @fn qcd_errorRoot(int test,int no,const char * const name,const char * const format,...)
 * @brief Evaluates a test on the root process and aborts if test is not 0
 * @details Checks whether "test"=0 on the root process and, if not, aborts the
 *     program gracefully with error number "no" after printing the "name"
 *     of the calling program and an error message to the error file from
 *     process 0. The message is formed using the "format" string and any
 *     additional arguments, exactly as in a printf() statement.
 * @param[in] test integer test condition
 * @param[in] no integer error number
 * @param[in] name 0-terminated string, that should indicate in what file/function
 *            the call was made
 * @param[in] format 0-terminated string with further error output, including 
 *            formatting as in printf()
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
void qcd_errorRoot(int test,int no,const char * const name,const char * const format,...);

/**
 * @fn qcd_error(int test,int no,const char * const name,const char * const format,...)
 * @brief Evaluates a test on all processes and aborts if test is not 0 everywhere.
 * @details Checks whether "test"=0 on every process and, if not, aborts the
 *     program gracefully with error number "no" after printing the "name"
 *     of the calling program and an error message to the error file from
 *     process 0. The message is formed using the "format" string and any
 *     additional arguments, exactly as in a printf() statement.
 *     This function must be called by all processes.
 * @param[in] test integer test condition
 * @param[in] no integer error number
 * @param[in] name 0-terminated string, that should indicate in what file/function
 *            the call was made
 * @param[in] format 0-terminated string with further error output, including 
 *            formatting as in printf()
 * @author Martin Luescher, minor adjustments Tomasz Korzec
 */
void qcd_error(int test,int no,const char * const name,const char * const format,...);

#endif
