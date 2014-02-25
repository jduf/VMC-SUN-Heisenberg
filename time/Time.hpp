#ifndef DEF_TIME
#define DEF_TIME

#include <ctime> 

class Time{
	public:
		Time():
			rawtime_(time(0)),
			time_(localtime(&rawtime_)) { }

		~Time(){};

		int day() const { return time_->tm_mday;}
		int month() const { return time_->tm_mon+1;}
		int year() const { return time_->tm_year+1900;}
		int hour() const { return time_->tm_hour;}
		int min() const { return time_->tm_min;}
		int sec() const { return time_->tm_sec;}

		bool limit_reached(time_t limit) const {
			if(time(0)-rawtime_ > limit){return true;}
			else{return false;}
		}

		void reset(){
			rawtime_ = time(0);
			time_ = localtime(&rawtime_);
		}

		time_t elapsed() const { return time(0)-rawtime_ ; }

	private:
		time_t rawtime_;
		struct tm* time_;
};
#endif
