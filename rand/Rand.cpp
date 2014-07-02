#include "Rand.hpp"

Rand::Rand(int len, int ij, int kl):
	c_ ((float)(  362436.0 / 16777216.0)),
	cd_((float)( 7654321.0 / 16777216.0)),
	cm_((float)(16777213.0 / 16777216.0)),
	i97_(96),
	j97_(32),
	len_(len),
	rvec_(new float[len_]),
	pos_(len_-1)
{
	init(ij,kl);
}

Rand::Rand(int len, Rand& seed):
	c_ ((float)(  362436.0 / 16777216.0)),
	cd_((float)( 7654321.0 / 16777216.0)),
	cm_((float)(16777213.0 / 16777216.0)),
	i97_(96),
	j97_(32),
	len_(len),
	rvec_(new float[len_]),
	pos_(len_-1)
{
	init(seed.get(31328),seed.get(30081));
}

Rand::Rand(int len):
	c_ ((float)(  362436.0 / 16777216.0)),
	cd_((float)( 7654321.0 / 16777216.0)),
	cm_((float)(16777213.0 / 16777216.0)),
	i97_(96),
	j97_(32),
	len_(len),
	rvec_(new float[len_]),
	pos_(len_-1)
{
	int seed1,seed2;
	sow(seed1,seed2);
#pragma omp critical
	{ init(seed1,seed2); }
}

Rand::Rand(int len, int thread):
	c_ ((float)(  362436.0 / 16777216.0)),
	cd_((float)( 7654321.0 / 16777216.0)),
	cm_((float)(16777213.0 / 16777216.0)),
	i97_(96),
	j97_(32),
	len_(len),
	rvec_(new float[len_]),
	pos_(len_-1)
{
	int seed1,seed2;
	sow(seed1,seed2);

	seed1+=thread;
	seed2+=thread;
#pragma omp critical
	{ init(seed1,seed2); }
}

void Rand::init(int ij, int kl){
	float s, t;
	int i, j, k, l, m;
	int ii, jj;

	if (ij < 0 || ij > 31328){ stop("ranmar: 1st seed must be between 0 and 31328\n");}
	if (kl < 0 || kl > 30081){ stop("ranmar: 2nd seed must be between 0 and 30081\n");}
	if (len_ < 0){ stop("ranmar: len < 0"); }

	i = (int)fmod(ij/177.0, 177.0) + 2;
	j = (int)fmod((double)ij, 177.0) + 2;
	k = (int)fmod(kl/169.0, 178.0) + 1;
	l = (int)fmod((double)kl, 169.0);

	for ( ii=0; ii<=96; ii++ ) {
		s = (float)0.0;
		t = (float)0.5;
		for ( jj=0; jj<=23; jj++ ) {
			m = (int)fmod( fmod((double)i*j,179.0)*k , 179.0 );
			i = j;
			j = k;
			k = m;
			l = (int)fmod( 53.0*l+1.0 , 169.0 );
			if ( fmod((double)l*m,64.0) >= 32){ s = s + t; }
			t = (float)(0.5 * t);
		}
		u_[ii] = s;
	}
}

Rand::~Rand() {
#pragma omp critical
	if(rvec_){ delete[] rvec_;}
}

void Rand::ranvec() {
	float uni;
	for (int ivec(0);ivec<len_;ivec++) {
		uni = u_[i97_] - u_[j97_];
		if ( uni < 0.0F ){ uni +=  1.0; }
		u_[i97_] = uni;
		i97_--;
		if ( i97_ < 0 ){ i97_ = 96; }
		j97_--;
		if ( j97_ < 0 ){ j97_ = 96; }
		c_ -=  cd_;
		if ( c_ < 0.0F ){ c_ = c_ + cm_; }
		uni -=  c_;
		if ( uni < 0.0F ){ uni = uni + 1.0; }
		rvec_[ivec] = uni;
	}
	pos_ = 0;
}

void Rand::stop(std::string msg) {
#pragma omp critical
	{ std::cerr << msg << std::endl; }
	exit(0);
}

void Rand::sow(int& seed1, int& seed2) {
	struct tm *tm_now;
	float s_sig, s_insig, maxs_sig, maxs_insig;
	time_t secs_now;
	int s, m, h, d, s1, s2;

	time(&secs_now);
	tm_now = localtime(&secs_now);

	s = tm_now->tm_sec + 1;
	m = tm_now->tm_min + 1;
	h = tm_now->tm_hour + 1;
	d = tm_now->tm_yday + 1;

	maxs_sig   = (float)(60.0 + 60.0/60.0 + 24.0/60.0/60.0 + 366.0/24.0/60.0/60.0);
	maxs_insig = (float)(60.0 + 60.0*60.0 + 24.0*60.0*60.0 + 366.0*24.0*60.0*60.0);

	s_sig      = (float)(s + m/60.0 + h/60.0/60.0 + d/24.0/60.0/60.0);
	s_insig    = (float)(s + m*60.0 + h*60.0*60.0 + d*24.0*60.0*60.0);

	s1 = (int)(s_sig   / maxs_sig   * 31328.0);
	s2 = (int)(s_insig / maxs_insig * 30081.0);

	seed1 = (s1+getpid())%31328;
	seed2 = s2;
}
