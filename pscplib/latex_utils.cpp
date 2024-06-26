/************************************************************************
* "latex_utils.c"
* Set of routines used for Latex tables
*
* JLP 
* Version 05/05/2010
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>    // exit()
#include <string.h>    // strcpy()
#include <ctype.h>                   /* isprint... */
#include "jlp_string.h" // jlp_trim_string (in jlplib/jlp_fits 

/* The prototypes of routines included here
*/
#include "latex_utils.h"

/*
#define DEBUG 
#define DEBUG_1 
*/
/**************************************************************************
* Read integer value in column #icol from b_data string
*
**************************************************************************/
int latex_read_ivalue(char *b_data, int *value, int icol)
{
int ival, status;
char buff[64];

*value = 0.;
status = latex_read_svalue(b_data, buff, icol);
if((status == 0) || (status == 3)) {
// Handle "\nnodata" value:
   if(status == 3) {
     *value = NO_DATA;
     status == 3;
     } else {
     ival = sscanf(buff, "%d", value);
/*
printf("latex_read_ivalue/buff=>%s< value=%d ival=%d\n", buff, *value, ival);
*/
     if(ival <= 0) status = 1;
     }
  }

return(status);
}
/**************************************************************************
* Read float value in column #icol from b_data string
*
* INPUT:
* iverbose : (i = 1 verbose if error)
*           (i > 1 verbose if error and even if no error)
**************************************************************************/
int latex_read_fvalue(char *b_data, float *value, int icol, int iverbose)
{
int ival, status;
char buff[64];

*value = 0.;
status = latex_read_svalue(b_data, buff, icol);
if((status == 0) || (status == 3)) {
// Handle "\nnodata" value:
   if(status == 3) {
     *value = NO_DATA;
     status == 3;
     } else {
     ival = sscanf(buff, "%f", value);
     if(ival <= 0) {
       if(iverbose > 1) printf("latex_read_fvalue/buff=>%s< value=%.2f ival=%d\n", buff, *value, ival);
       status = 1;
       }
     }
   }
return(status);
}
/**************************************************************************
* Read double value in column #icol from b_data string
*
* INPUT:
* iverbose : (i = 1 verbose if error)
*           (i > 1 verbose if error and even if no error)
**************************************************************************/
int latex_read_dvalue(char *b_data, double *value, int icol, int iverbose)
{
int ival, status;
char buff[64];

*value = 0.;
status = latex_read_svalue(b_data, buff, icol);
if((status == 0) || (status == 3)) {
// Handle "\nnodata" value:
   if(status == 3) {
     *value = NO_DATA;
     status == 3;
     } else {
     ival = sscanf(buff, "%lf", value);
     if(ival <= 0) {
        if(iverbose > 1) printf("latex_read_dvalue/buff=>%s< value=%.2f ival=%d\n", buff, *value, ival);
      status = 1;
      }
    }
  }
return(status);
}
/**************************************************************************
* Read string value in column #icol from b_data string
*
**************************************************************************/
int latex_read_svalue(char *b_data, char *value, int icol)
{
int ic, status, column_is_found;
char buff[NMAX], data[NMAX], *pc, nnodata_str[64];

strcpy(nnodata_str, "nodata");

*value = '\0';

strcpy(data, b_data);

pc = data;
data[NMAX-1] = '\0';
column_is_found = 0;
ic = 1;
buff[0] = '\0';
while(*pc && strncmp(pc,"\\\\",2)) {
  if(ic == icol) {
    column_is_found = 1;
    strcpy(buff,pc);
    break;
    }
  if(*pc == '&') {
    ic++;
    }
  pc++;
  }
*pc = '\0';
/* Return if column not found, or empty */
if(!buff[0]) return(-1);

/* Otherwise go on analysis: */
status = 1;
buff[NMAX-1] = '\0';
pc = buff;
while(*pc) {
  if(*pc == '&' || !strncmp(pc,"\\\\",2)) {
    *pc = '\0';
    strcpy(value,buff);
    if(*value) status = 0;
    break;
    }
  pc++;
  }

/* Removes '\r' (Carriage Return) if present: */
if(status == 0) {
pc = value;
while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
}

// Look for "\nnodata" string:
if(strstr(value, nnodata_str) != NULL) status = 3;

return(status);
}
/**************************************************************************
* Write a double value in column #icol to b_out string
*
* INPUT:
* b_data: current value of the line
* value: value to be written in b_data
* nber_of_decimal: number of decimals for the output format
*
* OUTPUT:
* b_out: updated value of the line b_data with the input value in Col.#icol
*
**************************************************************************/
int latex_write_dvalue(char *b_data, char *b_out, double value, int icol,
                       int nber_of_decimals)
{
int ic, column_is_found, istart, iend;
char data[360], *pc;
int i;

strcpy(data, b_data);

pc = data;
column_is_found = 0;
ic = 1;
i = 0;
while(*pc) {
  if(ic == icol && !column_is_found) {
    column_is_found = 1;
    istart = i;
    }
  else if(ic == icol+1) {
    iend = i-1;
    break;
    }
  if(*pc == '&') {
    ic++;
    }
  pc++;
  i++;
  }
/* Return if column not found, or empty */
if(istart == 0 || iend == istart) return(-1);

strcpy(b_out, b_data);

switch(nber_of_decimals) {
  case 1:
    sprintf(&b_out[istart],"%8.1f ",value);
    break;
  case 2:
    sprintf(&b_out[istart],"%8.2f ",value);
    break;
  case 3:
  default:
    sprintf(&b_out[istart],"%8.3f ",value);
    break;
  case 4:
    sprintf(&b_out[istart],"%9.4f ",value);
    break;
  }
strcpy(&b_out[istart+9],&b_data[iend]);

/*
printf("update_value/from >%s< to >%s< (value=%.2f)\n", b_data, b_out, value);
*/

return(0);
}
/************************************************************
* Get character string in the icol th column of a LaTeX table
*
* INPUT:
*  in_line: string corresponding to a full line of a Latex table
*  icol: column number
*
* OUTPUT:
*  item: character string corresponding to the column #icol
* 
************************************************************/
int latex_get_column_item(char *in_line, char *item, int icol, 
                          int verbose_if_error)
{
int ic;
char *pc, buffer[512];
 pc = in_line;

/* Get to the icol th column: */
 ic = 1;
 while(ic != icol) {
   while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3) 
        && strncmp(pc, "\\\\", 2)) pc++; 
    if(*pc == '&') {
      pc++;
      ic++;
     } else {
     if(verbose_if_error != 0) {
printf("verbose_if_error=%d\n", verbose_if_error);
        fprintf(stderr, "latex_get_column_item/Error: column #%d not found in >%s<\n",
             icol, in_line);
        }
     return(-1);
     }
 }

/* icol has been reached: */
 strcpy(buffer, pc);
 pc = buffer;
/* Go to the next column or the end of line: */
 while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3)
        && strncmp(pc, "\\\\", 2)) pc++; 
 *pc = '\0';
 strcpy(item, buffer);

return(0);
}
/************************************************************
* Set character string in the icol th column of a LaTeX table
*
* INPUT:
*  new_item: character string corresponding to the column #icol
*  new_item_len : length of the new item
*  icol: column number
*
* INPUT/OUTPUT:
*  in_line: string corresponding to a full line of a Latex table
* 
************************************************************/
int latex_set_column_item(char *in_line, int in_line_length,
                          char *new_item, int new_item_len, 
                          int icol, int verbose_if_error)
{
int ic, i1, i2, istart, i;
char *pc, buffer[512], out_line[512];;

// printf("set_column_item/icol=%d new_item_len=%d\n", icol, new_item_len);

 strcpy(out_line, in_line);
 pc = in_line;

/* Get to the icol th column and determine i1 the index of the first "&"
* of the icol th item: */
 i1 = 0;
 ic = 1;
 while(ic != icol) {
   while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3) 
        && strncmp(pc, "\\\\", 2)) { pc++; i1++; }
    if(*pc == '&') {
      pc++; i1++; 
      ic++;
     } else {
     if(verbose_if_error != 12340) {
printf("verbose_if_error=%d\n", verbose_if_error);
        fprintf(stderr, "latex_set_column_item/Error: column #%d not found in >%s<\n",
             icol, in_line);
         exit(-1);
        }
     return(-1);
     }
 }

// icol has been reached: go to next "&" 
/* Determine i1 the index of the last "&" of the icol th item: */
i2 = i1 + 1;
pc++;
/* Go to the next column or the end of line: */
 while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3)
        && strncmp(pc, "\\\\", 2)) {pc++; i2++;} 

/* Copy the new item to the icol th item: */
 for(i = 0; i < new_item_len; i++) {
    if(new_item[i] != '\0') {
       out_line[i1 + 1 + i] = new_item[i];
    } else {
     break;
    }
  }
/* DEBUG
printf("DEBUG/in_line=%s\n", in_line);
printf("DEBUG/out_line=%s\n, new_item=%s (len=%d)  i1=%d i2=%d *i2=%c\n", 
      out_line, new_item, new_item_len, i1, i2, in_line[i2]);
*/
/* Copy the end of the input line: */
istart = i1 + 1 + i;
 for(i = 0; i < 512; i++) {
    out_line[istart + i] = in_line[i2 + i];
    if(in_line[i2 + i] == '\0') break;
    }
/* DEBUG
printf("DEBUG/out_line=%s\n", out_line);
*/

// Copy for returning the argument:
 strncpy(in_line, out_line, in_line_length);

return(0);
}
/***********************************************************************
*
* Remove the i_col th column
*
* INPUT:
*  ilen : length of in_line (=input line)
*  i_col : column number (from 1...)
***********************************************************************/
int latex_remove_column(char *in_line, int i_col, int ilen)
{
char *pc, buffer1[256], buffer2[256];
int icol;

if(ilen > 256) {
  fprintf(stderr, "Fatal error: ilen=%d > 256! \n", ilen);
  exit(-1);
 }

// Exit if comments only:
if(in_line[0] == '\%') return(0);

strcpy(buffer1, in_line);
buffer1[255] = '\0';

/* Go to the i_col th column: */
pc = buffer1;
icol = 1;
while(*pc) {
while(*pc && (*pc != '&') &&
      strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
if(*pc == '&') {
  icol++;
// If column number is OK, stop scanning buffer1 here and copy remaining
// of the line to buffer2 (starting with "&")
  if(icol == i_col) {
   pc++;
   sprintf(buffer2, "%s", pc );
   buffer2[255] = '\0';
// Cut buffer1 after "&" location:
   pc--;
   *pc = '\0';
// Goto to next "&" in buffer2:
   pc = buffer2;
   while(*pc && (*pc != '&') &&
      strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
   sprintf(buffer2, "%s", pc );
   break;
   }
// End of line:
  } else {
  *pc = '\0';
  }
  pc++;
}

/* DEBUG:
printf("FIRST/remove_column: icol=%d\n in: %s \n buffer1: >%s< \n buffer2: >%s<\n", 
        i_col, in_line, buffer1, buffer2);
*/

// Cut last column of buffer1:
if(i_col == 1) {
strcpy(buffer1, "");
} else {
icol = 1;
pc = buffer1;
while(*pc) {
 while(*pc && (*pc != '&') &&
      strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
 if(*pc == '&') {
   icol++;
// If column number stop buffer1 here and copy remaining of the line to buffer2
   if(icol == i_col) {
// Cut after the "&" in buffer1:
    pc++;
    *pc = '\0';
    break;
    }
// End of line:
   } else {
   *pc = '\0';
   }
 pc++;
 }
}

/* DEBUG:
printf("SECOND/remove_column: buffer1: >%s< buffer2: >%s< \n out: %s %s\n",
        buffer1, buffer2, buffer1, buffer2);
*/

sprintf(in_line, "%s %s", buffer1, buffer2);

return(0);
}
/*************************************************************************
* Remove dollar in input string
*************************************************************************/
int latex_remove_dollar(char *in_string, char *out_string, int len_string)
{
int status = 0, i, j;

j = 0;
for(i = 0; i < len_string; i++) {
  if(in_string[i] == '\0')
    break;
// Remove dollar and also blank space...
  else if((in_string[i] != '$') && (in_string[i] != ' '))
    out_string[j++] = in_string[i];
  }
out_string[j] = '\0';

#ifdef DEBUG
// printf("remove_dollar/ in: >%s< out= >%s<\n", in_string, out_string);
#endif
return(status);
}
/*************************************************************************
* Add empty columns to obtain ncols_max 
* INPUT:
*  in_string
* OUTPUT:
*  out_string
*************************************************************************/
int latex_add_emptycols_to_ncols(char *in_string, char *out_string,
                                 int len_string_max, int ncols_max)
{
char *pc;
int ii, ncols;

// Copy in_string to out_string:
strcpy(out_string, in_string);
pc = out_string;
ncols = 1;
ii = 0;

// Scan the line and count the columns:
// skip \nodata with second condition..
while( (*pc != '\0') && ((*pc != '\\') || (*(pc+1) != '\\')) ) {
  if(*pc == '&') ncols++;
  pc++; ii++;
  }

// Add empty columns:
while(ncols < ncols_max) {
  *pc = ' '; pc++; ii++;
  *pc = '&'; pc++; ii++;
  *pc = ' '; pc++; ii++;
  ncols++;
  }

// Write "\\" at the end of the line
*pc = '\\'; pc++; ii++;
*pc = '\\'; pc++; ii++;
*pc = '\0';

// Check if ii has gone too far:
if(ii > len_string_max) {
   fprintf(stderr, "latex_add_emptycols_to_ncols/Fatal error ii=%d > len_string_max=%d \n",
           ii, len_string_max);
   exit(-1);
   }

// printf(" UUU3: ncols=%d >%s< ii=%d Len_max=%d\n", ncols, out_string,
//          ii, len_string_max);
return(0);
}

