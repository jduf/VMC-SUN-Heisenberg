#ifndef DEF_TIME
#define DEF_TIME

#include <ctime>
#include <chrono>
#include <string>

class Time{
	public:
		/*!Constructor*/
		Time(){ set(); }
		/*!Default destructor*/
		~Time() = default;
		/*{Forbidden*/
		Time(Time const&) = delete;
		Time(Time&&) = delete;
		Time& operator=(Time) = delete;
		/*}*/

		/*!Set to present time*/
		void set(){
			t0_ = std::chrono::steady_clock::now();
			tt_ = std::chrono::system_clock::to_time_t ( std::chrono::system_clock::now() );
			time_ = localtime(&tt_);
		}

		/*!Returns the current day*/
		int day() const { return time_->tm_mday; }
		/*!Returns the current month*/
		int month() const { return time_->tm_mon+1; }
		/*!Returns the current year*/
		int year() const { return time_->tm_year+1900; }
		/*!Returns the current hour*/
		int hour() const { return time_->tm_hour; }
		/*!Returns the current min*/
		int min() const { return time_->tm_min; }
		/*!Returns the current sec*/
		int sec() const { return time_->tm_sec; }
		/*!Returns the date*/
		std::string date(std::string s){
			char tmp[20];
			//std::strftime(tmp,20,"%G-%m-%d_%H:%M:%S",localtime(&rawtime_));
			std::string format("%F_%H"+s+"%M"+s+"%S");
			std::strftime(tmp,20,format.c_str(),localtime(&tt_));
			return tmp;
		}

		/*!Returns true if time limit (in second) has been reached*/
		bool limit_reached(time_t const& limit) const
		{ return time(0)>limit+tt_; }

		/*!Returns the elapsed from the instantiation or last call of set*/
		double elapsed() const { return (std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now() - t0_)).count(); }

	private:
		std::chrono::steady_clock::time_point t0_;
		time_t tt_;      //!< return value of time(0)
		struct tm* time_;//!< return value of localtime(time(0))
};
#endif
