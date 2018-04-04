/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/Timer.C,v $
 * $Revision: 1.4 $
 * $Author: vanvoor4 $
 * $Date: 2008-04-30 14:05:15 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.3  2007/08/21 15:56:32  vanvoor4
 * Added a missing '*'
 *
 * Revision 1.2  2007/01/24 14:49:58  vanvoor4
 * Fixed some bugs and added support for the 3 timers instead of just one.
 *
 * Revision 1.2  2006/11/16 20:37:45  vanvoor4
 * Profile timer uses SIGPROF not SIGALRM
 *
 * Revision 1.1  2006/04/11 20:30:32  vanvoor4
 * Initial checkin.
 *
 */

#include <iostream>
#include <cstring>
#include <errno.h>
#include <Timer.H>
using namespace ASCbase;

long Timer::A_timer_interval_sec = 0; 
long Timer::A_real_elapsed_sec = 0; 
long Timer::A_virt_elapsed_sec = 0;
long Timer::A_prof_elapsed_sec = 0;

Timer::Timer(int interval_sec_in)
{  
  // Currently the timers expire about 10 ms after reaching zero.
  // Probably because of only interrupting for timing 100 times a second.
  // In our case, 10ms is not that important if the timer_interval is 
  // reasonably large.
  if(interval_sec_in > 0) A_timer_interval_sec = interval_sec_in;
  else std::cerr << "The timer interval is less than 1 second; Skipping . .\n";

  A_fail = false;
  A_started = false;
}

bool 
Timer::start()
{
  struct itimerval tval;
  tval.it_value.tv_sec = A_timer_interval_sec;
  tval.it_value.tv_usec = 0;
  tval.it_interval.tv_sec = A_timer_interval_sec;
  tval.it_interval.tv_usec = 0;

  struct sigaction alrm_act;
  alrm_act.sa_sigaction = sig_handler;
  alrm_act.sa_flags = SA_SIGINFO | SA_RESTART;
  sigemptyset(&(alrm_act.sa_mask));

  if(sigaction(SIGALRM, &alrm_act, NULL)){
    std::cerr << "Unable to create SIGALRM for the profile timer\n";
    A_fail = true;
  }else if(sigaction(SIGVTALRM, &alrm_act, NULL)){
    std::cerr << "Unable to create SIGVTALRM for the profile timer\n";
    A_fail = true;
  }else if(sigaction(SIGPROF, &alrm_act, NULL)){
    std::cerr << "Unable to create SIGPROF for the profile timer\n";
    A_fail = true;
  }else if(setitimer(ITIMER_REAL, &tval, 0) == -1){
    std::cerr << "Unable to set the real timer: " << strerror(errno) << "\n";
    A_fail = true;
  }else if(setitimer(ITIMER_VIRTUAL, &tval, 0) == -1){
    std::cerr << "Unable to set the virtual timer: " << strerror(errno) << "\n";
    A_fail = true;
  }else if(setitimer(ITIMER_PROF, &tval, 0) == -1){
    std::cerr << "Unable to set the profile timer: " << strerror(errno) << "\n";
    A_fail = true;
  }

  if(!A_fail) A_started = true;
  return A_started;
}
 
bool 
Timer::get_itimer(double *rv, int which, long elapsed_sec)
{
  struct itimerval tval;
  if( getitimer(which, &tval) == -1){
    std::cerr << "Unable to get the current value of the timer: " 
              << strerror(errno) << "\n";
    *rv = -1.0;
    return false;
  }
  // timers count down from starting value
  *rv = elapsed_sec + A_timer_interval_sec - tval.it_value.tv_sec;
  *rv -= tval.it_value.tv_usec / 1000000.0;
  if(*rv < 0) *rv = 0;
  return true;
}

bool
Timer::get(double *real, double *virt, double *prof)
{
  if(A_started && !A_fail){
    get_itimer(real, ITIMER_REAL, A_real_elapsed_sec);
    get_itimer(virt, ITIMER_VIRTUAL, A_virt_elapsed_sec);
    get_itimer(prof, ITIMER_PROF, A_prof_elapsed_sec);
    return true;
  }
  return false;
}

void 
Timer::sig_handler(int signum, siginfo_t *siginfo, void *ucontext)
{
 if(signum == SIGALRM) A_real_elapsed_sec += A_timer_interval_sec;
 else if(signum == SIGVTALRM) A_virt_elapsed_sec += A_timer_interval_sec;
 else if(signum == SIGPROF) A_prof_elapsed_sec += A_timer_interval_sec;
}

std::string 
Timer::get_local_time(char *buf, size_t buf_len)
{
  struct timeval tv;
  struct timezone tz;
  if(gettimeofday(&tv, &tz) == -1){
    int errsv = errno;
    std::string tmp = "Unable get time using the \"gettimeofday\" function:  ";
    tmp += strerror(errsv);
    std::cerr << tmp << "\n";
    *buf = 0;
    return buf;
  }
  struct tm ltime;
  if(!localtime_r(&(tv.tv_sec), &ltime)){
    int errsv = errno;
    std::string tmp = "Unable to get the local time using \"localtime_r\":  ";
    tmp += strerror(errsv);
    std::cerr << tmp << "\n";
    *buf = 0;
    return buf;
  }
  if(my_strftime(buf, buf_len, "%c", &ltime) == 0) buf[79] = 0;
  return buf;
}
