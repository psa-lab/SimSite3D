/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/Timer.H,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-13 15:41:43 $
 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2008/04/30 14:05:15  vanvoor4
 * The timer class now uses sigaction in the place of signal.
 *
 * Revision 1.2  2007/01/24 14:49:31  vanvoor4
 * Fixed some long standing problems including checking for incorrect
 * signals and support for the 3 different timers.
 *
 * Revision 1.1  2006/04/11 20:30:36  vanvoor4
 * Initial checkin.
 *
 */

#ifndef _TIMER_HEADER_INCLUDED
#define _TIMER_HEADER_INCLUDED

#include <sys/time.h>
#include <signal.h>
#include <string>

namespace SimSite3D{
  //! A Linux profile timer for sequential programs
  /*! @ingroup common 
   * This class has been updated to use sigaction() in the place of signal(). 
   */
  class Timer{
  public:
    //! Constructor, sets the timer interval
    /*! Sets the timer interval and initial timer to timer_interval_in
     * 
     * @param interval_sec_in Value for timer interval an initial timer 
     */
    Timer(int interval_sec_in=10);
    
    //! Do nothing destructor
    ~Timer() {;}
    
    //! Starts profile timer
    /*!
     * @return True if timer was started, else false
     */
    bool start();

    //! Gets the elapsed time in seconds
    /*!
     * @param real Ptr to loc to store the elapsed real time
     * @param virt Ptr to loc to store the elapsed virtual time
     * @param prof Ptr to loc to store the elapsed profile time
     * @return True if !fail() and started()
     */
    bool get(double *real, double *virt, double *prof);
    
    //! Get the local time in a string format.
    /*! Does not necessarily belong in the Timer class, but currently this is
     * as good a place as any.
     * 
     * @param buf A Cstring used in the call to strftime
     * @param buf_len Length of the Cstring buf
     * @return A string holding the broken down format of the local time.
     */
    static std::string get_local_time(char *buf, size_t buf_len);

    //! True if an error had occured
    bool fail() { return A_fail; }
    
    //! True if timers have been started (i.e. start() has been called 
    //! successfully)
    bool started() { return A_started; }

    //! The sig handler used when the timer expires
    /*! 
     * @param signum The signal number passed to the sig handler
     * @param siginfo Pointer to a siginfo_t struct -- ignored by this function
     * @param ucontext Void pointer to a ucontext_t -- ignored by this function
     */
    static void sig_handler(int signum, siginfo_t *siginfo, void *ucontext);

  private:
    //! Privatized copy constructor
    /*! There is no valid reason to copy a timer class
     */
    Timer(const Timer &) {}

    //! Get the elapsed time for "which" timer in seconds
    bool get_itimer(double *rv, int which, long elapsed_sec);

    //! Wrapper to eliminate buggy gcc warning about using "%c" in strftime.
    /*! From the man page for strftime: 
     * Some buggy versions of gcc complain about the use of %c: warning: `%c'
     * yields only last 2 digits of year in some locales.  Of course
     * programmers are encouraged to use %c, it gives the preferred date
     * and time representation. One meets all kinds of strange obfuscations to
     * circumvent this gcc problem. A relatively clean one is to add this
     * intermediate function.
     *
     * @param s Character array of size max
     * @param max Size of the array s
     * @param fmt Format string
     * @param tm The tm structure to convert to a string.
     * @return The number of chars written to s if the time fits, else 0.
     */
    static size_t my_strftime(char *s, size_t max, const char *fmt,
	                      const struct tm *tm)
    { return strftime(s, max, fmt, tm); }
    
    bool A_fail;        //!< False implies no outstanding errors
    bool A_started;     //!< True implies timers are running

    static long A_timer_interval_sec; //!< Number of sec used to init a timer
    static long A_real_elapsed_sec; //!< real elapsed time not including timer
    static long A_virt_elapsed_sec; //!< virtual elapsed time not including timer
    static long A_prof_elapsed_sec; //!< profile elapsed time not including timer
  };
}
#endif
