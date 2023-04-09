/********************************************************************
* jlp_calib_table.h
*
* JLP
* Version 03/03/2023
*********************************************************************/
#ifndef __jlp_calib_table_h // BOF sentry
#define  __jlp_calib_table_h

class JLP_CalibTab {

public:
    JLP_CalibTab(char *calib_table_fname0); 
    ~JLP_CalibTab() {}; 

// Look for next line with double star measurement:
    int NextLineWithMeasurements(char *full_line, 
                                 double *WDS_alpha, double *WDS_delta);

private:
   char calib1_fname[128];
   FILE *fp_calib1;
   int iline1;

};

#endif //EOF sentry
