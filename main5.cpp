#include <iostream>
#include <vector>
#include "dt_heating_model.hpp"
#include <cmath> 

/**
Определяет время перегрева в 75*С при различных обратных тяговых токах
*/

template <class T > 
T mabs( T data )
{
    return data < 0 ? -data : data;
}

std::vector<double> curents  {500,1000,1500};//,2000,2500,3000,3500,4000,4500,5000,5500,6000};

unsigned heats( double current, double max_temp, Sci::dt_heat_transfer& dt_temp_calc )
{
    unsigned heating_sec=0;
    for ( ; heating_sec < 200001 and dt_temp_calc.oil_temp_degC() < max_temp; ++heating_sec  )

        dt_temp_calc.calc( current, 1 );

    return heating_sec;
}



int main ()
{

    std::cout << "curretn,06_light_time,06_time,04_time " << std::endl; 	

    for ( auto a : curents )
    {
        Sci::dt dt_06_1000_light( "DT-0.6-1000",
			40, 390, 0.0172, 288, 10.472, 0.58, 
			70, 480, 0.019, 
			24.3, 1670, 
			47.0, 480.0, 0.76, 0.258, 1.3, 0.86, 0.2585, 0.8, 2000 );

	Sci::dt dt_06_1000( "ДТ-0.6-1000",
		    69, 390, 0.0172, 288, 10.872, 0.58, 
		    120, 480, 0.283, 
		    41.3, 1670, 
		    81.0, 480.0, 0.76, 0.3, 1.3, 0.86, 0.258, 0.9, 
		    2000 );

	Sci::dt dt_04_1500( "ДТ-0.4-1500",
		    104, 390, 0.0172, 288, 10.872, 0.767, 
		    182, 480, 0.374, 
		    63.3, 1670, 
		    122, 480.0, 1.005, 0.396, 1.719, 1.137, 0.341, 0.9, 
		    3000 );

	Sci::heat_transfer_data htd( 100.0, 100.0, 100.0, 5.35, 8.7, 800.0, 273.0, 1, 0);
	Sci::dt_temp_value hc( 273.0, 273.0, 273.0, 273.0 );

        Sci::dt_heat_transfer dt_06_temp_calc_light ( dt_06_1000_light, htd, hc );
        Sci::dt_heat_transfer dt_06_temp_calc ( dt_06_1000, htd, hc );
        Sci::dt_heat_transfer dt_04_temp_calc ( dt_04_1500, htd, hc );

	
	std::cout <<  a << ", " << heats ( a, 75, dt_06_temp_calc_light )<< "," << heats ( a, 75, dt_06_temp_calc )<< "," << heats ( a, 75, dt_04_temp_calc ) << std::endl; 	

    }

    return 0;
}

