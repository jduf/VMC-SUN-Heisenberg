#ifndef DEF_TIME
#define DEF_TIME

#include <ctime> 

class Time{
	public:
		Time(){
			time_t t(time(0));
			now=localtime(&t);
		};
		~Time(){};

		int day() const { return now->tm_mday;}
		int month() const { return now->tm_mon+1;}
		int year() const { return now->tm_year+1900;}
		int hour() const { return now->tm_hour;}
		int min() const { return now->tm_min;}
		int sec() const { return now->tm_sec;}

	private:
		struct tm* now;
};
#endif
