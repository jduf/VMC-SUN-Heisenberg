#ifndef DEF_TIME
#define DEF_TIME

#include <ctime> 

class Time{
	public:
		/*!Constructor*/
		Time(){set();}
		/*!Destructor*/
		~Time(){}
		/*!Set to present time*/
		void set(){
			rawtime_ = time(0);
			time_ = localtime(&rawtime_);
		}

		/*!Returns the current day*/
		int day() const { return time_->tm_mday;}
		/*!Returns the current month*/
		int month() const { return time_->tm_mon+1;}
		/*!Returns the current year*/
		int year() const { return time_->tm_year+1900;}
		/*!Returns the current hour*/
		int hour() const { return time_->tm_hour;}
		/*!Returns the current min*/
		int min() const { return time_->tm_min;}
		/*!Returns the current sec*/
		int sec() const { return time_->tm_sec;}

		/*!Returns true if limit*/
		bool limit_reached(time_t limit) const {
			return (elapsed() > limit)?true:false;
		}

		/*!Returns the elapsed from the instantiation or last call of set*/
		time_t elapsed() const { return time(0)-rawtime_; }

	private:
		/*!Forbids copy*/
		Time(Time const& l);
		/*!Forbids assignment*/
		Time& operator=(Time const& l);

		time_t rawtime_;	//!< return value of time(0)
		struct tm* time_;	//!< return value of localtime(time(0))
};
#endif
