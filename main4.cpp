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

void calc_trains_packet( double current, std::vector< std::pair<double, double>  > & value, int max_trains_packet = 20, double move_time = 60.0, double free_time = 60.0 )
{
    Sci::dt dt_06_1000(40, 390, 0.0172, 288, 7.472, 0.58, 70, 480, 0.019, 24.3, 1670, 47.0, 480.0, 0.76, 0.258, 1.3, 0.86, 0.2585, 0.8, 2000 );
    Sci::heat_transfer_data htd( 310.0, 310.0, 310.0, 8.35, 8.7, 800.0, 273.0);
    Sci::heat_condition hc( 273.0, 273.0, 273.0, 273.0 );
    Sci::dt_heat_transfer dt_temp_calc_(dt_06_1000, htd, hc );

    int calc_step_in_sec = 20;

    int move_time_iter = move_time * 60 / calc_step_in_sec;
    int free_time_iter = free_time * 60 / calc_step_in_sec;

std::cout << "move_time_iter: " << move_time_iter << std::endl;
std::cout << "free_time_iter: " << free_time_iter << std::endl;

    for ( int step = 0; step <= max_trains_packet; ++step )	
    {
        //heating
	for ( int i = 0; i <= move_time_iter; ++i )
	{
	    dt_temp_calc_.calc( current, calc_step_in_sec );
	    value.push_back( std::make_pair( dt_temp_calc_.oil_temp_K() - htd.external_temp(), calc_fail_rate_correction ( current/dt_06_1000.nominal_current(), dt_temp_calc_.oil_temp() ) ) );
	}
	    //cooling
	for ( int i = 0; i <= free_time_iter; ++i )
	{
	    dt_temp_calc_.calc( 0, calc_step_in_sec );
	    value.push_back( std::make_pair( dt_temp_calc_.oil_temp_K() - htd.external_temp(), calc_fail_rate_correction ( current/dt_06_1000.nominal_current(), dt_temp_calc_.oil_temp() ) ) );
	}
    }
}



int main ()
{
	std::vector< std::pair< double, double > > free60, free50, free40, free30, free20;
	calc_trains_packet( 1535, free60, 10, 60, 60 );
	calc_trains_packet( 1535, free50, 10, 60, 50 );
	calc_trains_packet( 1535, free40, 10, 60, 40 );
	calc_trains_packet( 1535, free30, 10, 60, 30 );
	calc_trains_packet( 1535, free20, 10, 60, 20 );


    std::cout << "Toil_free60; Alfa_free60; Toil_free50; Alfa_free50; Toil_free40; Alfa_free40; Toil_free30; Alfa_free30; Toil_free20; Alfa_free20" << std::endl;
    for (int i = 0; i < free60.size(); ++i )
    {
        std::cout << free60[i].first << "; " << free60[i].second << "; " << free50[i].first << "; " << free50[i].second << "; " << free40[i].first << "; " << free40[i].second << "; ";
	std::cout << free30[i].first << "; " << free30[i].second << "; " << free20[i].first << "; " << free20[i].second << std::endl;

    }













return 0;
}

