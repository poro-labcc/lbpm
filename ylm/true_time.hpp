//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
// Classe True_Time
//
//   Function that retirns the current time.
//
//________________________________________________________
//A.Z. - 02/13 => Creation
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef TRUE_TIME_HPP
#define TRUE_TIME_HPP

#include <sys/time.h>


class True_Time{
   public:
   
     True_Time(void){};
     
     double inline operator () (void){
        gettimeofday(&_t, NULL);
        return ( _t.tv_sec + _t.tv_usec/1000000.0 );
     }


   private:
      struct timeval _t;
};


#endif // TRUE_TIME_HPP
