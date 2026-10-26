/** @file qcd_io.h
 *  @brief header file for qcd_io.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 *  @todo: rename get* --> read*
 */

#ifndef H_QCD_IO
#define H_QCD_IO 1 

/* constants */
/**
 * @def qcd_GF_LIME
 * @brief Identifier for gauge fields in ILDG format
 * @todo: format description
 */
#define qcd_GF_LIME 0

/**
 * @def qcd_GF_OPENQCD
 * @brief Identifier for gauge fields in openQCD (export) format
 * @todo format description
 */
#define qcd_GF_OPENQCD 1

/**
 * @def qcd_GF_SFMPHMC
 * @brief Identifier for gauge fields in Hubert's lexicographical format (little endian)
 * @details Used e.g. in Marina's SFMPHMC code or in the sfqcd branch of openQCD
 * @todo format description
 */
#define qcd_GF_SFMPHMC 2
/**
 * @def qcd_GF_SFMPHMCb
 * @brief Identifier for gauge fields in Hubert's lexicographical format (big endian)
 * @details Used e.g. in Marina's SFMPHMC code or in the sfqcd branch of openQCD
 * @todo format description
 */
#define qcd_GF_SFMPHMCb 3

/**
 * @def qcd_PROP_LIME
 * @brief Identifier for propagators in ILDG format
 * @todo format description
 */
#define qcd_PROP_LIME 0

/**
 * @def qcd_PROP_OPENQCD
 * @brief Identifier for propagators in openQCD format
 * @todo format description
 */
#define qcd_PROP_OPENQCD 1

/**
 * @def qcd_PROP_SFMPHMC
 * @brief Identifier for propagators in SFMPHMC format (little endian)
 * @todo format description
 */
#define qcd_PROP_SFMPHMC 2

/**
 * @def qcd_PROP_SFMPHMCb
 * @brief Identifier for propagators in SFMPHMC format (big endian)
 * @todo format description
 */
#define qcd_PROP_SFMPHMCb 3

/**
 * @def qcd_PROP_HMC
 * @brief Identifier for propagators in the format of ETMC's inverter
 * @todo format description
 */
#define qcd_PROP_HMC 4

/**
 * @def qcd_PROP_HMCV
 * @brief Identifier for propagators in the format of ETMC's inverter
 * @details The 12 spinor fields contained in a propagator are stored in 
 *          separate files.
 * @todo format description
 */
#define qcd_PROP_HMCV 5

/**
 * @def qcd_PROP_CMI
 * @brief Identifier for propagators in Chris Michael's format
 * @todo format description
 */
#define qcd_PROP_CMI 6

/**
 * @def qcd_SOURCE_LIME
 * @brief Identifier for spinor fields in IDLG format
 * @todo format description
 */
#define qcd_SOURCE_LIME 103


/* prototypes */

/**
 * @fn char* qcd_readFile(char* fname, long *len)
 * @brief Reads an ASCII file into a character array.
 * @details The file with name given by fname is read and the resulting 
 *          character array returned. The Length of the file is written into
 *          the long integer pointed to by len.
 * @param[in] fname contains the file name
 * @param len *len will be overwritten with the length of the file in bytes
 * @return a pointer to an array of characters.
 */
char* qcd_readFile(char* fname, long *len);

/**
 * @fn int qcd_GFformatID(char *format)
 * @brief Given the gauge field format as a string, the corresponding integer is returned.
 *        If format does not represent a known  gauge field format, -1 is returned.
 * @param[in] format a 0-terminated string describing the gauge field format, e.g. 
 *            "qcd_GF_LIME"
 * @return the integer that is used internally to represent this format
 */
int qcd_GFformatID(char *format);

/**
 * @fn int qcd_getLimeMessage(char *fname, qcd_geometry *geo, char **message_ptr)
 * @brief Reads the "lime message" from an ILDG gauge configuration file.
 * @details Process 0 reads the lime message from file fname using the lime 
 *          library. The message
 *          is broadcasted to all other processes and a pointer to it is written
 *          to *message_ptr
 * @param fname 0-terminated string, the file-name of the ILDG gauge configuration
 * @param geo pointer to an initialized geometry structure
 * @param message_ptr address of a pointer to the read message
 * @return 0 on success, 1 on failure
 * @todo geo is not really needed (numprocs and myid can be obtained from MPI)
 */
int qcd_getLimeMessage(char *fname, qcd_geometry *geo, char **message_ptr);

/**
 * @fn int qcd_getGaugeLime(char *fname, qcd_gaugeField *u)
 * @brief Reads a gauge configuration in ILDG format.
 * @details The gauge field in the file with the name given by fname is read
 *          using parallel I/O and stored in the gauge field u.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 * @todo output arguments should come first
 */
int qcd_getGaugeLime(char *fname, qcd_gaugeField *u);

/**
 * @fn int qcd_getGaugeOpenQCD(char *fname, qcd_gaugeField *u)
 * @brief Reads a gauge configuration in openQCD format.
 * @details The gauge field in the file with the name given by fname is read
 *          using parallel I/O and stored in the gauge field u.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 * @todo output arguments should come first
 */
int qcd_getGaugeOpenQCD(char *fname, qcd_gaugeField *u);

/**
 * @fn int qcd_getGaugeOpenQCDSeq(char *fname, qcd_gaugeField *u)
 * @brief Reads a gauge configuration in openQCD format.
 * @details The gauge field in the file with the name given by fname is read
 *          using sequential I/O and stored in the gauge field u.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 * @todo output arguments should come first
 */
int qcd_getGaugeOpenQCDSeq(char *fname, qcd_gaugeField *u);

/**
 * @fn int qcd_getGaugeSFMPHMC(char *fname, qcd_gaugeField *u, int littleEndian)
 * @brief Reads a gauge configuration in SFMPHMC format.
 * @details The gauge field in the file with the name given by fname is read
 *          using parallel I/O and stored in the gauge field u.
 *          If littleEndian is 1 it is assumed that the file is in little endian,
 *          if it is 0, it is assumed that the file is in big endian.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @param[in] littleEndian specifies the byte ordering of the file
 * @return 0 on success, 1 on failure
 * @todo output arguments should come first
 */
int qcd_getGaugeSFMPHMC(char *fname, qcd_gaugeField *u, int littleEndian);

/**
 * @fn int qcd_getGaugeField(char *fname, int type, qcd_gaugeField *u)
 * @brief Reads a gauge configuration.
 * @details The gauge field in the file with the name given by fname, which has
 *          the format given by type is read
 *          using parallel I/O and stored in the gauge field u.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in] type identifier of the file format. E.g. qcd_GF_LIME
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 * @todo output arguments should come first
 */
int qcd_getGaugeField(char *fname, int type, qcd_gaugeField *u);

/**
 * @fn int qcd_writeGaugeLime(char *fname, qcd_gaugeField *u, char *message)
 * @brief Writes a gauge configuration in ILDG format.
 * @details The gauge field u is written to a file  with the name given by fname
 *          using parallel I/O. The meta-information that is stored in the 
 *          lime message must be provided in the string message.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] message 0-terminated string containing the meta-data
 * @return 0 on success, 1 on failure
 */
int qcd_writeGaugeLime(char *fname, qcd_gaugeField *u, char *message);

/**
 * @fn int qcd_writeGaugeOpenQCD(char *fname, qcd_gaugeField *u)
 * @brief Writes a gauge configuration in openQCD's (export) format.
 * @details The gauge field u is written to a file  with the name given by fname
 *          using parallel I/O. 
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_writeGaugeOpenQCD(char *fname, qcd_gaugeField *u);

/**
 * @fn int qcd_writeGaugeSFMPHMC(char *fname, qcd_gaugeField *u, int littleEndian)
 * @brief Writes a gauge configuration in SFMPHMC format.
 * @details The gauge field u is written to a file  with the name given by fname
 *          using parallel I/O. If littleEndian is 1, the file is written 
 *          with little endian byte ordering, if it is 0, big endian ordering is used.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] littleEndian specifies the byte ordering to use
 * @return 0 on success, 1 on failure
 */
int qcd_writeGaugeSFMPHMC(char *fname, qcd_gaugeField *u, int littleEndian);

/**
 * @fn int qcd_writeGaugeField(char *fname, int type, qcd_gaugeField *u, ...)
 * @brief Writes a gauge configuration in a specified format.
 * @details The gauge field u is written to a file with the name given by fname
 *          using parallel I/O. The file format is specified by type, e.g. 
 *          qcd_GF_LIME. In case of the ILDG format, the lime-message has to 
 *          be provided as an optional parameter.
 *          All processes must call this function simultaneously.
 * @param[in] fname 0-terminated string providing the file-name
 * @param[in] type one of the supported format identifiers
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] ... in case of the qcd_GF_LIME format, a 0-terminated string
 *                containing the lime-message has to be provided in addition
 * @return 0 on success, 1 on failure
 */
int qcd_writeGaugeField(char *fname, int type, qcd_gaugeField *u, ...);

/**
 * @fn int qcd_getPropagatorCMI(char *fname, qcd_propagator *p)
 * @brief Reads a propagator from a CMI file.
 * @details The propagator in the file with the name given by fname, which has
 *          the CMI format, is read
 *          using parallel I/O and stored in the propagator p.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getPropagatorCMI(char *fname, qcd_propagator *p);

/**
 * @fn int qcd_getPropagatorHMC(char *fname, qcd_propagator *p)
 * @brief Reads a propagator from an HMC file.
 * @details The propagator in the file with the name given by fname, which has
 *          the HMC format, is read
 *          using parallel I/O and stored in the propagator p.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getPropagatorHMC(char *fname, qcd_propagator *p);

/**
 * @fn int qcd_getPropagatorLime(char *fname, qcd_propagator *p)
 * @brief Reads a propagator from an ILDG file.
 * @details The propagator in the file with the name given by fname, which has
 *          the ILDG format, is read
 *          using parallel I/O and stored in the propagator p.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @author Giannis Koutsou
 * @todo output parameters should come first
 */
int qcd_getPropagatorLime(char *fname, qcd_propagator *p);

/**
 * @fn int qcd_getPropagator(char *fname, int type, qcd_propagator *p)
 * @brief Reads a propagator from a file with a specified format.
 * @details The propagator in the file with the name given by fname, which has
 *          a format given by type (e.g. qcd_PROP_LIME), is read
 *          using parallel I/O and stored in the propagator p.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] type the file-format identifier
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getPropagator(char *fname, int type, qcd_propagator *p);

/**
 * @fn int qcd_getVectorCMI(char *fname, int nu, int c2, qcd_vector *v)
 * @brief Reads a spinor field from a CMI propagator file.
 * @details One of the 12 spinor fields contained in a propagator is read from 
 *          the CMI-file with name given by fname. The result is stored in v.
 *          if the propagator in the file has components p->D[x][mu][nu][c1][c2],
 *          then the "column" \n
 *          v->D[x][mu][c1] = p->D[x][mu][nu][c1][c2], for all x,mu,c1\n 
 *          is read.
 *          Parallel I/O is used, so
 *          all processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] nu Dirac component of propagator which to read
 * @param[in] c2 color component of propagator which to read
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getVectorCMI(char *fname, int nu, int c2, qcd_vector *v);

/**
 * @fn int qcd_getVectorHMC(char *fname, int nu, int c2, qcd_vector *v)
 * @brief Reads a spinor field from a HMC propagator file.
 * @details One of the 12 spinor fields contained in a propagator is read from 
 *          the HMC-file with name given by fname. The result is stored in v.
 *          if the propagator in the file has components p->D[x][mu][nu][c1][c2],
 *          then the "column" \n
 *          v->D[x][mu][c1] = p->D[x][mu][nu][c1][c2], for all x,mu,c1\n 
 *          is read.
 *          Parallel I/O is used, so
 *          all processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] nu Dirac component of propagator which to read
 * @param[in] c2 color component of propagator which to read
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getVectorHMC(char *fname, int nu, int c2, qcd_vector *v);

/**
 * @fn int qcd_getVectorHMCV(char *fname, int nu, int c2, qcd_vector *v)
 * @brief Reads a spinor field from HMCV propagator files.
 * @details One of the 12 spinor fields contained in a propagator is read from 
 *          the HMCV-file with name given by fname. The result is stored in v.
 *          if the propagator in the file has components p->D[x][mu][nu][c1][c2],
 *          then the "column" \n
 *          v->D[x][mu][c1] = p->D[x][mu][nu][c1][c2], for all x,mu,c1\n 
 *          is read.
 *          Parallel I/O is used, so
 *          all processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] nu Dirac component of propagator which to read
 * @param[in] c2 color component of propagator which to read
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @todo output parameters should come first
 */
int qcd_getVectorHMCV(char *fname, int nu, int c2, qcd_vector *v);

/**
 * @fn int qcd_getVectorLime(char *fname, qcd_vector *v)
 * @brief Reads a spinor field from a spinor field file in ILDG format.
 * @details The spinor field in the file with the name given by fname, which has
 *          the ILDG format, is read
 *          using parallel I/O and stored in the spinor field v.
 *          All processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_getVectorLime(char *fname, qcd_vector *v);

/**
 * @fn int qcd_getVectorOpenQCD(char *fname, qcd_vector *v)
 * @brief Reads a spinor field from a spinor field file in openQCD format.
 * @details The spinor field in the file with the name given by fname, which has
 *          the openQCD format, is read
 *          using parallel I/O and stored in the spinor field v.
 *          All processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_getVectorOpenQCD(char *fname, qcd_vector *v);

/**
 * @fn int qcd_getVectorSFMPHMC(char *fname, qcd_vector *v, int getLE)
 * @brief Reads a spinor field from a spinor field file in SFMPHMC format.
 * @details The spinor field in the file with the name given by fname, which has
 *          the SFMPHMC format, is read
 *          using parallel I/O and stored in the spinor field v. If getLE==1,
 *          the file is assumed to be stored in little endian byte ordering, if 
 *          getLE==0, big endian ordering is assumed.
 *          All processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in] getLE specifies the byte ordering
 * @return 0 on success, 1 on failure
 */
int qcd_getVectorSFMPHMC(char *fname, qcd_vector *v, int getLE);

/**
 * @fn int qcd_getVector(char *fname, int type, int nu, int c2, qcd_vector *v)
 * @brief Reads a spinor field from a file in a specified format.
 * @details The spinor field in the propagator file with the name given by fname, which has
 *          a format given by type (e.g. #qcd_PROP_LIME), is read
 *          using parallel I/O and stored in the spinor field v. If the file 
 *          contains the spinor p, the component\n
 *          v->D[x][mu][c1] = p->D[x][mu][nu][c1][c2] \n
 *          are read, for a fixed, given nu and c1.
 *          All processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in] type file format specifier
 * @param[in] nu dirac component of propagator to read
 * @param[in] c2 color component of propagator to read
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_getVector(char *fname, int type, int nu, int c2, qcd_vector *v);

/**
 * @fn int qcd_writePropagatorCMI(char *fname, qcd_propagator *p)
 * @brief Writes a propagator to a CMI file.
 * @details The propagator p is written to the file with the name given 
 *          by fname in the #qcd_PROP_CMI format, using parallel I/O.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo if read* parameter order changes, this too...
 */
int qcd_writePropagatorCMI(char *fname, qcd_propagator *p);

/**
 * @fn qcd_writePropagatorHMCV(char *fname, qcd_propagator *p)
 * @brief Writes a propagator to a HMCV file.
 * @details The propagator p is written to the file with the name given 
 *          by fname in the #qcd_PROP_HMCV format, using parallel I/O.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo if read* parameter order changes, this too...
 */
int qcd_writePropagatorHMCV(char *fname, qcd_propagator *p);

/**
 * @fn int qcd_writePropagator(char *fname, int type, qcd_propagator *p)
 * @brief Writes a propagator to a file with a specified format.
 * @details The propagator p is written to the file with the name given by 
 *          fname, in a format given by type (e.g. #qcd_PROP_LIME),
 *          using parallel I/O.
 *          All processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] type the file-format identifier
 * @param[in] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 * @todo if read* parameter order changes, this too...
 */
int qcd_writePropagator(char *fname, int type, qcd_propagator *p);

/**
 * @fn int qcd_writeVectorCMI(char *fname, qcd_vector *v)
 * @brief Writes a spinor field to a spinor file.
 * @details The spinor field v is written to 
 *          the file with name given by fname in the #qcd_PROP_CMI format.
 *          Parallel I/O is used, so
 *          all processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @todo if read* parameter order changes, this too
 */
int qcd_writeVectorCMI(char *fname, qcd_vector *v);

/**
 * @fn int qcd_writeVectorLime(char *fname, int type, qcd_vector *v)
 * @brief Writes a spinor field to a file in ILDG format.
 * @details The spinor field v is written to  the file with the name given 
 *          by fname, which has
 *          the ILDG format (either #qcd_PROP_LIME or #qcd_SOURCE_LIME).
 *          Parallel I/O is used, so all processes must call this 
 *          function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in] type file format specifier
 * @param[in] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @author Giannis Koutsou
 */
int qcd_writeVectorLime(char *fname, int type, qcd_vector *v);

/**
 * @fn int qcd_writeVectorOpenQCD(char *fname, qcd_vector *v)
 * @brief Writes a spinor field to a spinor field file in openQCD format.
 * @details The spinor field v is written to the file with the name given 
 *          by fname, which has the qcd_PROP_OPENQCD format.
 *          Parallel I/O is used, so 
 *          All processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_writeVectorOpenQCD(char *fname, qcd_vector *v);

/**
 * @fn int qcd_writeVector(char *fname, int type, qcd_vector *v)
 * @brief Writes a spinor field to a file with specified format.
 * @details The spinor field v is written to a file with
 *          with name given by fname, in a format specified by type
 *          (e.g. #qcd_PROP_LIME). 
 *          Parallel I/O is used, so
 *          all processes must call this function simultaneously.
 * @param[in] fname file-name given as a 0-terminated string 
 * @param[in] type file format specifier
 * @param[in] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 * @todo if read* parameter order changes, this too
 */
int qcd_writeVector(char *fname, int type, qcd_vector *v);

/**
 * @fn int qcd_writeVectorSFMPHMC(char *fname, qcd_vector *v, int writeLE);
 * @brief Writes a spinor field to a file in SFMPHMC format.
 * @details The spinor field v is written to the file with the name 
 *          given by fname, which has
 *          the #qcd_PROP_SFMPHMC or #qcd_PROP_SFMPHMCb format, depending on 
 *          the byte ordering flag writeLE (1=little endian, 0=big endian).
 *          Parallel I/O is used, 
 *          all processes must call this function simultaneously.
 * @param[in] fname file name as a 0-terminated string
 * @param[in] v pointer to an initialized qcd_vector structure
 * @param[in] writeLE specifies the byte ordering
 * @return 0 on success, 1 on failure
 */
int qcd_writeVectorSFMPHMC(char *fname, qcd_vector *v, int writeLE);

/**
 * @fn int qcd_getSpinorComponents3d(char *fname, int nmodes, qcd_spinorComponent3d *evec)
 * @brief Reads 3d spinor-component fields from a file.
 * @brief nmodes time-slice spinor-component fields are read from file *fname
 *        and stored in evec[0]...evec[nmodes-1].
 *        All processes of a 3d communicator must call this function simultaneously.
 *        If more than 1 3d communicator is involved, the 3d spinors of the 
 *        different communicators will read the same content, unless their 
 *        fname is different, which is allowed.
 * @param[in] fname file name given as a 0-terminated string.
 * @param[in] nmodes number of spinor-component fields to read
 * @param[in,out] evec pointer to an initialized qcd_spinorComponent3d structure
 * @return  0 on success, 1 on failure
 * @todo parameter order
 */
int qcd_getSpinorComponents3d(char *fname, int nmodes, qcd_spinorComponent3d *evec);

/**
 * @fn int qcd_writeSpinorComponents3d(char *fname, int nmodes, qcd_spinorComponent3d *evec)
 * @brief Writes 3d spinor-component fields to a file.
 * @brief nmodes time-slice spinor-component fields evec[0]...evec[nmodes-1] are 
 *        written to a file *fname.
 *        All processes of a 3d communicator must call this function simultaneously.
 *        If more than one 3d communicator is involved, the different communicators
 *        should have different fname, otherwise the behaviour is undefined.
 * @param[in] fname file name given as a 0-terminated string.
 * @param[in] nmodes number of spinor-component fields to read
 * @param[in] evec pointer to an initialized qcd_spinorComponent3d structure
 * @return  0 on success, 1 on failure
 * @todo parameter order
 */
int qcd_writeSpinorComponents3d(char *fname, int nmodes, qcd_spinorComponent3d *evec);
#endif
