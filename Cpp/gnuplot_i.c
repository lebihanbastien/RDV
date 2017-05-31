/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.c
  @author   N. Devillard
  @date     Sep 1998
  @version  $Revision: 2.10 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/
/*--------------------------------------------------------------------------*/

/*
    $Id: gnuplot_i.c,v 2.10 2003/01/27 08:58:04 ndevilla Exp $
    $Author: ndevilla $
    $Date: 2003/01/27 08:58:04 $
    $Revision: 2.10 $
 */

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/
#include "gnuplot_i.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#ifdef WIN32
#include <io.h>
#endif // #ifdef WIN32

/*---------------------------------------------------------------------------
                                Defines
 ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
                          Prototype Functions
 ---------------------------------------------------------------------------*/

/**
 * Creates a temporary file name for writing
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 *
 * @return char const * Pointer to file name string.
 */
char const * gnuplot_tmpfile(gnuplot_ctrl * handle);

/**
 * Plot a temporary file.
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 * @param tmp_filename
 * @param title
 */
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle,
                           char const* tmp_filename,
                           char const* title,
                           char const* ls,
                           char const* lt,
                           char const* lw,
                           int lc);

void gnuplot_plot_atmpfile_3d(gnuplot_ctrl * handle,
                              char const* tmp_filename,
                              char const* title,
                              char const* ls,
                              char const* lt,
                              char const* lw,
                              int lc);

/*---------------------------------------------------------------------------
                            Function codes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
/*--------------------------------------------------------------------------*/

gnuplot_ctrl * gnuplot_init(void)
{
    gnuplot_ctrl *  handle ;
    int i;

#ifndef WIN32
    if (getenv("DISPLAY") == NULL) {
        fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
    }
#endif // #ifndef WIN32


    /*
     * Structure initialization:
     */
    handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;
    gnuplot_setcolor(handle, 0);


    handle->gnucmd = popen("gnuplot", "w") ;
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
        free(handle) ;
        return NULL ;
    }

    for (i=0;i<GP_MAX_TMP_FILES; i++)
    {
        handle->tmp_filename_tbl[i] = NULL;
    }
    return handle;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
/*--------------------------------------------------------------------------*/

void gnuplot_close(gnuplot_ctrl * handle)
{
    int     i ;

    if (pclose(handle->gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    if (handle->ntmp) {
        for (i=0 ; i<handle->ntmp ; i++) {
            remove(handle->tmp_filename_tbl[i]) ;
            free(handle->tmp_filename_tbl[i]);
            handle->tmp_filename_tbl[i] = NULL;

        }
    }
    free(handle) ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Sends a command to an active gnuplot session.
  @param    handle Gnuplot session control handle
  @param    cmd    Command to send, same as a printf statement.

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  @code
  gnuplot_cmd(g, "plot %d*x", 23.0);
  gnuplot_cmd(g, "plot %.18e * cos(%.18e * x)", 32.0, -3.0);
  @endcode

  Since the communication to the gnuplot process is run through
  a standard Unix pipe, it is only unidirectional. This means that
  it is not possible for this interface to query an error status
  back from gnuplot.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...)
{
    va_list ap ;

    va_start(ap, cmd);
    vfprintf(handle->gnucmd, cmd, ap);
    va_end(ap);

    fputs("\n", handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Change the plotting style of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_style Plotting-style to use (character string)
  @return   void

  The provided plotting style is a character string. It must be one of
  the following:

  - lines
  - points
  - linespoints
  - impulses
  - dots
  - steps
  - errorbars
  - boxes
  - boxeserrorbars
 */
/*--------------------------------------------------------------------------*/

void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style)
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        strcpy(h->pstyle, "points") ;
    } else {
        strcpy(h->pstyle, plot_style) ;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Change the plotting color of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_color Plotting-color to use (int)
  @return   void
 */
/*--------------------------------------------------------------------------*/

void gnuplot_setcolor(gnuplot_ctrl * h, int plot_style)
{
    switch(plot_style)
    {
        case 1:
        strcpy(h->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(h->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(h->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(h->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(h->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(h->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(h->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(h->pcolor, "black") ;
        break;
        case 9:
        strcpy(h->pcolor, "gray50") ;
        break;
        default:
        strcpy(h->pcolor, "default") ;
        break;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_xlabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set xlabel \"%s\"", label) ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set ylabel \"%s\"", label) ;
}

void gnuplot_set_zlabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set zlabel \"%s\"", label) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_resetplot(gnuplot_ctrl * h)
{
    int     i ;
    if (h->ntmp) {
        for (i=0 ; i<h->ntmp ; i++) {
            remove(h->tmp_filename_tbl[i]) ;
            free(h->tmp_filename_tbl[i]);
            h->tmp_filename_tbl[i] = NULL;

        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Plots a 2d graph from a list of doubles.
  @param    handle  Gnuplot session control handle.
  @param    d       Array of doubles.
  @param    n       Number of values in the passed array.
  @param    title   Title of the plot.
  @return   void

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  @code
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot_x(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_x(
    gnuplot_ctrl    *   handle,
    double          *   d,
    int                 n,
    char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || d==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n ; i++) {
      fprintf(tmpfd, "%.18e\n", d[i]);
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title, "lines", "1", "2", 1);
    return ;
}



/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a 2d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

  @code
    gnuplot_ctrl    *h ;
    double          x[50] ;
    double          y[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        x[i] = (double)(i)/10.0 ;
        y[i] = x[i] * x[i] ;
    }
    gnuplot_plot_xy(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_xy(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    char const           *   title,
    char const      * ls,
    char const      * lt,
    char const      * lw,
    int lc)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title, ls, lt, lw, lc);
    return ;
}


void gnuplot_plot_xyz(
    gnuplot_ctrl    * handle,
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * title,
    char const      * ls,
    char const      * lt,
    char const      * lw,
    int lc)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e %.18e\n", x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile_3d(handle,tmpfname,title, ls, lt, lw, lc);
    return ;
}


void gnuplot_fplot_xy(
    double          * x,
    double          * y,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e \n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}


void gnuplot_fplot_xyz(
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e, %.18e, %.18e\n", x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}

void gnuplot_fplot_txyz(
    double          * t,
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e %.18e %.18e\n", t[i], x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}



/*-------------------------------------------------------------------------*/
/**
  @brief    Open a new session, plot a signal, close the session.
  @param    title   Plot title
  @param    style   Plot style
  @param    label_x Label for X
  @param    label_y Label for Y
  @param    x       Array of X coordinates
  @param    y       Array of Y coordinates (can be NULL)
  @param    n       Number of values in x and y.
  @return

  This function opens a new gnuplot session, plots the provided
  signal as an X or XY signal depending on a provided y, waits for
  a carriage return on stdin and closes the session.

  It is Ok to provide an empty title, empty style, or empty labels for
  X and Y. Defaults are provided in this case.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_once(
  char    *   title,
  char    *   style,
  char    *   label_x,
  char    *   label_y,
  double  *   x,
  double  *   y,
  int         n
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || n<1) return ;

  if ((handle = gnuplot_init()) == NULL) return ;
  if (style!=NULL) {
      gnuplot_setstyle(handle, style);
  } else {
      gnuplot_setstyle(handle, "lines");
  }
  if (label_x!=NULL) {
      gnuplot_set_xlabel(handle, label_x);
  } else {
      gnuplot_set_xlabel(handle, "X");
  }
  if (label_y!=NULL) {
      gnuplot_set_ylabel(handle, label_y);
  } else {
      gnuplot_set_ylabel(handle, "Y");
  }
  if (y==NULL) {
      gnuplot_plot_x(handle, x, n, title);
  } else {
      gnuplot_plot_xy(handle, x, y, n, title, "lines", "1", "2", 1);
  }
  printf("press ENTER to continue\n");
  while (getchar()!='\n') {}
  gnuplot_close(handle);
  return ;
}

void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(handle, "%s %.18e * x + %.18e title \"%s\" with %s",
                  cmd, a, b, title, handle->pstyle) ;

    handle->nplots++ ;
    return ;
}


void gnuplot_plot_equation(
    gnuplot_ctrl    *   h,
    char            *   equation,
    char            *   title
)
{
    char const *    cmd    = (h->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(h, "%s %s title \"%s\" with %s",
                  cmd, equation, title, h->pstyle) ;
    h->nplots++ ;
    return ;
}


int gnuplot_write_x_csv(
    char const * fileName,
    double const * d,
    int n,
    char const * title)
{
    int     i;
    FILE*   fileHandle;

    if (fileName==NULL || d==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e\n", i, d[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xy_csv(
    char const *        fileName,
    double const    *   x,
    double const    *   y,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || y==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%.18e, %.18e\n", x[i], y[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_multi_csv(
    char const *        fileName,
    double const    **  xListPtr,
    int                 n,
    int                 numColumns,
    char const      *   title)
{
    int     i;
    int     j;
    FILE*   fileHandle;

    if (fileName==NULL || xListPtr==NULL || (n<1) || numColumns <1)
    {
        return -1;
    }

    for (j=0;j<numColumns;j++)
    {
        if (xListPtr[j] == NULL)
        {
            return -1;
        }
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e", i, xListPtr[0][i]) ;
        for (j=1;j<numColumns;j++)
        {
            fprintf(fileHandle, ", %.18e", xListPtr[j][i]) ;
        }
        fprintf(fileHandle, "\n");
    }

    fclose(fileHandle) ;

    return 0;
}

char const * gnuplot_tmpfile(gnuplot_ctrl * handle)
{
    static char const * tmp_filename_template = "gnuplot_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    int                 tmp_filelen = strlen(tmp_filename_template);

#ifndef WIN32
    int                 unx_fd;
#endif // #ifndef WIN32

    assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);

    /* Open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return NULL;
    }

    tmp_filename = (char*) malloc(tmp_filelen+1);
    if (tmp_filename == NULL)
    {
        return NULL;
    }
    strcpy(tmp_filename, tmp_filename_template);

#ifdef WIN32
    if (_mktemp(tmp_filename) == NULL)
    {
        return NULL;
    }
#else // #ifdef WIN32
    unx_fd = mkstemp(tmp_filename);
    if (unx_fd == -1)
    {
        return NULL;
    }
//    close(unx_fd);

#endif // #ifdef WIN32

    handle->tmp_filename_tbl[handle->ntmp] = tmp_filename;
    handle->ntmp ++;
    return tmp_filename;
}

/*
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    if(strcmp(handle->pcolor,"default")==0) gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lw 2", cmd, tmp_filename,title, handle->pstyle);
    else gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lt rgb lw 2 \"%s\"", cmd, tmp_filename,title, handle->pstyle, handle->pcolor);

    handle->nplots++ ;
    return ;
}
*/


void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title, char const* ls, char const* lt, char const* lw, int lc)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    switch(lc)
    {
        case 1:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(handle->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(handle->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(handle->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(handle->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(handle->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(handle->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(handle->pcolor, "black") ;
        break;
        case 9:
        strcpy(handle->pcolor, "gray50") ;
        break;
        default:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
    }

    gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lt \"%s\" lw %s lc rgb \"%s\"", cmd, tmp_filename,title, ls, lt, lw, handle->pcolor);
    handle->nplots++ ;
    return ;
}



void gnuplot_plot_atmpfile_3d(gnuplot_ctrl * handle,
                              char const* tmp_filename,
                              char const* title,
                              char const* ls,
                              char const* lt,
                              char const* lw,
                              int lc)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "splot";
    title                  = (title == NULL)      ? "(none)" : title;

    switch(lc)
    {
        case 1:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(handle->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(handle->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(handle->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(handle->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(handle->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(handle->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(handle->pcolor, "black") ;
        break;
        case 9:
        strcpy(handle->pcolor, "gray50") ;
        break;
        default:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
    }
    gnuplot_cmd(handle, "set view 60,75");
    gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lt \"%s\" lw %s lc rgb \"%s\"", cmd, tmp_filename, title, ls, lt, lw, handle->pcolor);
    handle->nplots++ ;
    return ;
}

/* vim: set ts=4 et sw=4 tw=75 */
