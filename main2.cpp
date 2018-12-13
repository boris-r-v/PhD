#include <iostream>
#include <vector>
#include "sum_temp.hpp"
#include <cmath> 

std::vector<double> ef_curent_at_time_inerv  {1535,1426,1395,1294};


double calc_fail_rate_correction( float load_coef, float temp )
{
    float ret = 0;
    if (load_coef * temp >= 30)
    { /*Зона высоких нагрузок*/
	ret = 0.008434*std::exp(0.048281*temp)*std::exp(4.08637*load_coef);
    }
    else
    {	/*Зона низких нагрузок*/
	ret = -4.137+0.01*(temp+273)+2.42*load_coef;
    }
    return ret > 0 ? ret : 0.1;
}

void calc_trains_packet( double current, std::vector< std::pair<double, double>  > & value, int max_trains_packet = 20, double move_time = 60.0 )
{
    Sci::dt dt_06_1000(40, 390, 0.0172, 288, 7.472, 0.58, 70, 480, 0.019, 24.3, 1670, 47.0, 480.0, 0.76, 0.258, 1.3, 0.86, 0.2585, 0.8, 2000 );
    Sci::heat_transfer_data htd( 310.0, 310.0, 310.0, 8.35, 8.7, 800.0, 273.0);
    Sci::heat_condition hc( 273.0, 273.0, 273.0, 273.0 );
    Sci::dt_heat_transfer dt_temp_calc_(dt_06_1000, htd, hc );

    int calc_num_steps = 20.0;

    double calc_step = move_time / (double)calc_num_steps;

    for ( int step = 0; step <= max_trains_packet; ++step )	
    {
        //heating
	for ( int i = 0; i <= calc_num_steps; ++i )
	{
	    //heating
	    dt_temp_calc_.calc( current, calc_step );
	    value.push_back( std::make_pair( dt_temp_calc_.oil_temp_K() - htd.external_temp(), calc_fail_rate_correction ( current/dt_06_1000.nominal_current(), dt_temp_calc_.oil_temp() ) ) );
	    
	}
	    //cooling
	for ( int i = 0; i <= calc_num_steps; ++i )
	{
	    dt_temp_calc_.calc( 0, calc_step );
	    value.push_back( std::make_pair( dt_temp_calc_.oil_temp_K() - htd.external_temp(), calc_fail_rate_correction ( current/dt_06_1000.nominal_current(), dt_temp_calc_.oil_temp() ) ) );
	}
    }
}



int main ()
{

	std::vector<std::pair<double, double> > one, min6, min8, min10, min12;	
	calc_trains_packet( 1870, one, 30, 20 );
	calc_trains_packet( 1535, min6, 10, 60 );
	calc_trains_packet( 1429, min8, 10, 60  );
	calc_trains_packet( 1359, min10, 10, 60  );
	calc_trains_packet( 1294, min12, 10, 60 );
	
    std::cout << "Toil_alone; Alfa_alone; Toil_6min; Alfa_6min; Toil_8min; Alfa_8min; Toil_10min; Alfa_10min; Toil_12min; Alfa_12min" << std::endl;

    for (int i = 0; i < min6.size(); ++i )
        std::cout << one[i].first << "; " << one[i].second << "; " << min6[i].first << "; " << min6[i].second << "; " << min8[i].first << "; " << min8[i].second << "; " << min10[i].first << "; " << min10[i].second << "; " << min12[i].first << "; " << min12[i].second << std::endl;














return 0;
}

