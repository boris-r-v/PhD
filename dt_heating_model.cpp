#include "dt_heating_model.hpp"
#include <cmath>
static const double Csb (5.6e-8);

void Sci::heat_transfer_data::change_wind( )
{
    std::cout << "change_wind: old wind: " << wind_ << "m/s" << std::endl;
    wind_ = distribution_(generator_);
    std::cout << "change_wind: new wind: " << wind_ << "m/s" << std::endl;

}

double Sci::heat_transfer_data::get_Pr( double _degC )
{
    _degC -= 273;
    if ( _degC <= 0 ) return 866;
					//=0,0002291667*B14^4-0,0314166667*B14^3+1,7620833333*B14^2-52,9083333*B14+866
    if ( 0 < _degC and _degC <= 40 ) return (0.0002291667*std::pow(_degC, 4 ) - 0.0314166667*std::pow(_degC, 3 ) + 1.7620833333*std::pow(_degC, 2 ) - 52.9083333*_degC + 866 );
					//=2,5553613053*10^-6*G19^4 - 0,0010518065*G19^3 + 0,169326049*G19^2 - 13,1306351981*G19 + 459,987
    if ( 40 < _degC and _degC <=120 ) return (0.000001*2.5553613053*std::pow(_degC, 4 ) - 0.0010518065*std::pow(_degC, 3 ) + 0.169326049*std::pow(_degC, 2 ) - 13.1306351981*_degC + 459.99 );
    if ( _degC > 120 ) return 34.9;
}

double Sci::heat_transfer_data::get_epsilon( double _degK_fluid, double _degK_solid )
{
    return std::pow ( get_Pr(_degK_fluid) / get_Pr ( _degK_solid ), 0.25 );
}

    //h_body_air_=2.1+2.5*std::pow((Tbody - Tair),0.25)*std::pow((Tair/Tbody), 0.25);	-	это из дисера
    //h_body_air_ = 2.1+1.2834*std::log(Tbody-272.9)+1.51;		//это апрксимированная
    //std::cout << "Tbody: " << Tbody << std::endl;
void Sci::heat_transfer_data::h_parametrs_recalc( dt_temp_value const& dt  )
{
    h_body_air_=2.5*std::pow((dt.body_temp() - external_temp_),0.25)*std::pow((293.0/dt.body_temp()), 0.25) + 4.2*wind_;
    //double base (0.4019*dt.oil_temp() + 110.09);
    
    double base (0.5255*dt.oil_temp() - 75.898);

    h_coil_oil_ = base * get_epsilon( dt.oil_temp(), dt.coil_temp() );
    h_oil_body_ = base * get_epsilon( dt.oil_temp(), dt.body_temp() );
    h_oil_core_ = base * get_epsilon( dt.oil_temp(), dt.core_temp() );  

    //std::cout << "coil-oil: " << get_epsilon( dt.oil_temp(), dt.coil_temp() ) << std::endl;
    //std::cout << "oil-core: " << get_epsilon( dt.oil_temp(), dt.core_temp() ) << std::endl;
    //std::cout << "oil-body: " << get_epsilon( dt.oil_temp(), dt.body_temp() ) << std::endl;
}


Sci::dt::dt( 	std::string const& title,
		double coil_M, double coil_HC, double coil_R, double coil_CSA, double coil_L, double coil_S,
		double core_M, double core_HC, double core_S,
		double oil_M, double oil_HC,
		double body_M, double body_HC, double body_oil_S, double body_air_S, double body_radiation_S, double body_ground_S, double body_sun_S, double body_grayness,
		double  nominal_current ):
	title_( title ),
	coil_M_(coil_M), coil_HC_(coil_HC), coil_R_(coil_R), coil_CSA_(coil_CSA), coil_L_(coil_L), coil_S_(coil_S),
	core_M_(core_M), core_HC_(core_HC), core_S_(core_S),
	oil_M_(oil_M), oil_HC_(oil_HC),
	body_M_(body_M), body_HC_(body_HC), body_oil_S_(body_oil_S), body_air_S_(body_air_S), body_radiation_S_(body_radiation_S), body_ground_S_(body_ground_S), body_sun_S_(body_sun_S), body_grayness_(body_grayness),
	nominal_current_( nominal_current )
{
}


Sci::heat_transfer_data::heat_transfer_data(double h_coil_oil, double h_oil_body, double h_oil_core, double h_body_air, double h_body_ground, double sun_radiation, double external_temp, double cloud, double wind ):
    h_coil_oil_(h_coil_oil),  h_oil_body_(h_oil_body),  h_oil_core_(h_oil_core),  
    h_body_air_(h_body_air),  h_body_ground_(h_body_ground), sun_radiation_(sun_radiation),  
    external_temp_(external_temp),
    cloud_( 1 - cloud ),
    generator_(),
    distribution_( wind, wind/3 ),
    wind_ (  distribution_(generator_) )
{
}
Sci::dt_temp_value::dt_temp_value():
     oil_temp_( NAN ),  coil_temp_( NAN ),  core_temp_( NAN ),  body_temp_( NAN )
{
}


Sci::dt_temp_value::dt_temp_value( double oil_temp, double coil_temp, double core_temp, double body_temp ):
     oil_temp_(oil_temp),  coil_temp_(coil_temp),  core_temp_(core_temp),  body_temp_(body_temp)
{
}

bool Sci::dt_temp_value::valid() const
{
    return !isnan( oil_temp_ ) and !isnan( coil_temp_ ) and !isnan( core_temp_ ) and !isnan (body_temp_);
}

Sci::dt_heat_transfer::dt_heat_transfer( dt const& dt, heat_transfer_data& htd, dt_temp_value& hc ):
    dt_(dt),
    htd_(htd),
    hc_( hc )
{
}

void Sci::dt_heat_transfer::set_temp_value( dt_temp_value const& _dt )
{
    hc_ = _dt;
}

void Sci::dt_heat_transfer::calc( double _current, double _calculate_sec )
{
    /*Пересчитаем коэффициенты теплоотдачи на начальную температуру вычислений*/
    htd_.h_parametrs_recalc( hc_ );

    rk45_error err;
    current( _current );    
    double dt ( 0.1 );
    double steps_number( 1 + _calculate_sec/dt );
    double t0 ( hc_.oil_temp() );

//    std::cout  << "Calc Step: " << dt << "(sec)" << std::endl;
//    std::cout << "STart: " << *this << std::endl;
    for ( int i = 0; i <= steps_number; ++i )
    {
	rk45( dt, err );
//	std::cout << "STEP: " << i << " " << *this << " oil error: " << err.oil << std::endl;
    }


}

void Sci::dt_heat_transfer::rk45( double _dt, rk45_error& _error )
{
    dt_temp_value const& hc_const = hc_;    
    double kcoil0 = _dt * coil_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double koil0  = _dt * oil_ode(  hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double kcore0 = _dt * core_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double kbody0 = _dt * body_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );

    double kcoil1_c ( hc_const.coil_temp() + kcoil0/2 );
    double koil1_c  ( hc_const.oil_temp() + koil0/2 );
    double kcore1_c ( hc_const.core_temp() + kcore0/2 );
    double kbody1_c ( hc_const.body_temp() + kbody0/2 );
    double kcoil1 = _dt * coil_ode( kcoil1_c, koil1_c, kcore1_c, kbody1_c );
    double koil1  = _dt * oil_ode ( kcoil1_c, koil1_c, kcore1_c, kbody1_c );
    double kcore1 = _dt * core_ode( kcoil1_c, koil1_c, kcore1_c, kbody1_c );
    double kbody1 = _dt * body_ode( kcoil1_c, koil1_c, kcore1_c, kbody1_c );

    double kcoil2_c ( hc_const.coil_temp() + (kcoil0 + kcoil1)/4 );
    double koil2_c  ( hc_const.oil_temp() +  (koil0 + koil1)/4   );
    double kcore2_c ( hc_const.core_temp() + (kcore0 + kcore1)/4 );
    double kbody2_c ( hc_const.body_temp() + (kbody0 + kbody1)/4 );
    double kcoil2 = _dt * coil_ode( kcoil2_c, koil2_c, kcore2_c, kbody2_c );
    double koil2  = _dt * oil_ode ( kcoil2_c, koil2_c, kcore2_c, kbody2_c );
    double kcore2 = _dt * core_ode( kcoil2_c, koil2_c, kcore2_c, kbody2_c );
    double kbody2 = _dt * body_ode( kcoil2_c, koil2_c, kcore2_c, kbody2_c );

    double kcoil3_c ( hc_const.coil_temp() - kcoil1 + 2*kcoil2 );
    double koil3_c  ( hc_const.oil_temp()  - koil1 + 2*koil2 );
    double kcore3_c ( hc_const.core_temp() - kcore1 + 2*kcore2 );
    double kbody3_c ( hc_const.body_temp() - kbody1 + 2*kbody2);
    double kcoil3 = _dt * coil_ode( kcoil3_c, koil3_c, kcore3_c, kbody3_c );
    double koil3  = _dt * oil_ode ( kcoil3_c, koil3_c, kcore3_c, kbody3_c );
    double kcore3 = _dt * core_ode( kcoil3_c, koil3_c, kcore3_c, kbody3_c );
    double kbody3 = _dt * body_ode( kcoil3_c, koil3_c, kcore3_c, kbody3_c );

    double kcoil4_c ( hc_const.coil_temp() + (7*kcoil0 + 10*kcoil1 + kcoil3)/27 );
    double koil4_c  ( hc_const.oil_temp() + (7*koil0 + 10*koil1 + koil3)/27 );
    double kcore4_c ( hc_const.core_temp() + (7*kcore0 + 10*kcore1 + kcore3)/27 );
    double kbody4_c ( hc_const.body_temp() + (7*kbody0 + 10*kbody1 + kbody3)/27 );
    double kcoil4 = _dt * coil_ode( kcoil4_c, koil4_c, kcore4_c, kbody4_c );
    double koil4  = _dt * oil_ode ( kcoil4_c, koil4_c, kcore4_c, kbody4_c );
    double kcore4 = _dt * core_ode( kcoil4_c, koil4_c, kcore4_c, kbody4_c );
    double kbody4 = _dt * body_ode( kcoil4_c, koil4_c, kcore4_c, kbody4_c );

    double kcoil5_c ( hc_const.coil_temp() + (28*kcoil0 - 128*kcoil1 + 546*kcoil2 + 54*kcoil3 - 378*kcoil4)/625 );
    double koil5_c  ( hc_const.oil_temp() +  (28*koil0 - 128*koil1 + 546*koil2 + 54*koil3 - 378*koil4)/625 );
    double kcore5_c ( hc_const.core_temp() + (28*kcore0 - 128*kcore1 + 546*kcore2 + 54*kcore3 - 378*kcore4)/625 );
    double kbody5_c ( hc_const.body_temp() + (28*kbody0 - 128*kbody1 + 546*kbody2 + 54*kbody3 - 378*kbody4)/625 );
    double kcoil5 = _dt * coil_ode( kcoil5_c, koil5_c, kcore5_c, kbody5_c );
    double koil5  = _dt * oil_ode ( kcoil5_c, koil5_c, kcore5_c, kbody5_c ); 
    double kcore5 = _dt * core_ode( kcoil5_c, koil5_c, kcore5_c, kbody5_c );
    double kbody5 = _dt * body_ode( kcoil5_c, koil5_c, kcore5_c, kbody5_c );

    double coil_4ord ( hc_const.coil_temp() + ( kcoil0 + 4*kcoil2 + kcoil3 )/6 ); 
    double oil_4ord  (hc_const.oil_temp()  +  ( koil0 +  4*koil2 +  koil3  )/6 );  
    double core_4ord (hc_const.core_temp() +  ( kcore0 + 4*kcore2 + kcore3 )/6 );   
    double body_4ord (hc_const.body_temp() +  ( kbody0 + 4*kbody2 + kbody3 )/6 );   

    double coil_5ord ( hc_const.coil_temp() + ( 14*kcoil0 + 35*kcoil3 + 162*kcoil4 + 125*kcoil5 )/336 ); 
    double oil_5ord  (hc_const.oil_temp()  +  ( 14*koil0 +  35*koil3  + 162*koil4  + 125*kcoil5 )/336 );  
    double core_5ord (hc_const.core_temp() +  ( 14*kcore0 + 35*kcore3 + 162*kcore4 + 125*kcore5 )/336 );   
    double body_5ord (hc_const.body_temp() +  ( 14*kbody0 + 35*kbody3 + 162*kbody4 + 125*kbody5 )/336 );   

    
    hc_.coil_temp( coil_5ord );
    hc_.oil_temp ( oil_5ord  );
    hc_.core_temp( core_5ord );
    hc_.body_temp( body_5ord );

    _error.coil = coil_4ord - coil_5ord;
    _error.oil = oil_4ord - oil_5ord;
    _error.core = core_4ord - core_5ord;
    _error.body = body_4ord - body_5ord;

}
//(10.872*0.0172*Io^2*(1+0.004*(Tcoil-293))/288)

double Sci::dt_heat_transfer::coil_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    return ( ( current()*current()*(1+0.004*(_coil_temp - 293))*dt_.coil_L()*dt_.coil_R()/dt_.coil_CSA() ) - (htd_.h_coil_oil()*dt_.coil_S()*(_coil_temp - _oil_temp) ) ) / ( dt_.coil_M()*dt_.coil_HC() );
}
double Sci::dt_heat_transfer::oil_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    return ( (htd_.h_coil_oil()*dt_.coil_S()*(_coil_temp - _oil_temp) ) - (htd_.h_oil_body()*dt_.body_oil_S()*(_oil_temp - _body_temp)) - (htd_.h_oil_core()*dt_.core_S()*(_oil_temp - _core_temp))  
	    - (htd_.h_body_air()*0.1*(_oil_temp - htd_.external_temp() ) ) ) / ( dt_.oil_M()*dt_.oil_HC() );
}
double Sci::dt_heat_transfer::core_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    return ( htd_.h_oil_core()*dt_.core_S()*(_oil_temp - _core_temp) ) / ( dt_.core_M()*dt_.core_HC() );
}
double Sci::dt_heat_transfer::body_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    double colling ( htd_.h_body_air()*dt_.body_air_S()*(_body_temp - htd_.external_temp()) +	//ковекция с корпуса
		    htd_.h_body_ground()*dt_.body_ground_S()*(_body_temp - htd_.external_temp() ) + //кондукция с корпуса
		    dt_.body_grayness()*dt_.body_radiation_S()*Csb*(std::pow(_body_temp,4) - std::pow( htd_.external_temp(),4 ) ) //излучение с корпуса
		    );

    return ( (htd_.h_oil_body()*dt_.body_oil_S()*(_oil_temp - _body_temp)) + (dt_.body_grayness()*dt_.body_sun_S()*htd_.sun_radiation()*htd_.cloud()) - colling ) / ( dt_.body_M()*dt_.body_HC() );
}


std::ostream& Sci::operator << ( std::ostream& _s, dt const& _d )
{
    _s << "Тип ДТ: " << _d.title() << std::endl;
    _s << "Параметры корпуса: "<< std::endl;
    _s << " Общая площадь поверхности (излучения): " << _d.body_radiation_S() <<" (м^2)" << std::endl;
    _s << " Площадь поверхности контакта с воздухом: " << _d.body_air_S() <<" (м^2)" << std::endl;
    _s << " Площадь поверхности контакта с землей: " << _d.body_ground_S() <<" (м^2)" << std::endl;
    _s << " Площадь поверхности контакта с маслом: " << _d.body_oil_S() <<" (м^2)" << std::endl;
    _s << " Площадь поверхности нагрева солнечной радиацией: " << _d.body_sun_S() <<" (м^2)" << std::endl;
    _s << " Масса: " << _d.body_M() <<" (кг)" << std::endl;
    _s << " Теплоемкость: " << _d.body_HC() <<" (Дж/К)" << std::endl;

    _s << "Параметры масла: "<< std::endl;
    _s << " Площадь поверхности контакта с корпусом: " << _d.body_oil_S() <<" (м^2)" << std::endl;
    _s << " Площадь поверхности контакта с сердечником: " << _d.core_S() <<" (м^2)" << std::endl;
    _s << " Масса: " << _d.oil_M() <<" (кг)" << std::endl;
    _s << " Теплоемкость: " << _d.oil_HC() <<" (Дж/К)" << std::endl;

    _s << "Параметры сердечника: "<< std::endl;
    _s << " Площадь поверхности контакта с маслом: " << _d.core_S() <<" (м^2)" << std::endl;
    _s << " Mасса: " << _d.core_M() <<" (кг)" << std::endl;
    _s << " Теплоемкость: " << _d.core_HC() <<" (Дж/К)" << std::endl;

    _s << "Параметры обмотки: "<< std::endl;
    _s << " Площадь поверхности контакта с маслом: " << _d.coil_S() <<" (мм^2)" << std::endl;
    _s << " Удельное сопротивление материала: " << _d.coil_R() <<" (Ом*мм^2/м)" << std::endl;
    _s << " Площадь поперечного сечения: " << _d.coil_CSA() <<" (мм^2)" << std::endl;
    _s << " Общая длина: " << _d.coil_L() <<" (м)" << std::endl;
    _s << " Mасса: " << _d.coil_M() <<" (кг)" << std::endl;
    _s << " Теплоемкость: " << _d.coil_HC() <<" (Дж/К)" << std::endl;


    return _s;
}

std::ostream& Sci::operator << ( std::ostream& _s, heat_transfer_data const& _h )
{
    _s << std::endl;
    _s << "Параметры теплопередачи:" << std::endl;
    _s << "  Коэф.теплоотдачи обмотка-масло: " << _h.h_coil_oil() << " (Вт/(м^2*K)" << std::endl;
    _s << "  Коэф.теплоотдачи масло-корпус: " << _h.h_oil_body() << " (Вт/(м^2*K)" << std::endl;
    _s << "  Коэф.теплоотдачи масло-сердечник: " << _h.h_oil_core() << " (Вт/(м^2*K)" << std::endl;
    _s << "  Коэф.теплоотдачи корпус-воздух: " << _h.h_body_air() << " (Вт/(м^2*K)" << std::endl;
    _s << "  Коэф.теплоотдачи корпус-земля: " << _h.h_body_ground() << " (Вт/(м^2*K)" << std::endl;
    _s << "  Мошность солнечного излучения: " << _h.sun_radiation() << " (Вт)" << std::endl;
    _s << "  Снижение солнеч.изл. из-за облачности: " << _h.cloud() << std::endl;
    _s << "  Температурв охлаждающей среды: " << _h.external_temp() << " (*K) {" << _h.external_temp()-273 << "*C}" << std::endl;
    return _s;
}

std::ostream& Sci::operator << (std::ostream& _s, dt_heat_transfer const& _ht)
{
    heat_transfer_data const& htd = _ht.get_heat_transfer_data();
    _s << "{Tcoil="<< _ht.coil_temp_degC() << ", Toil=" << _ht.oil_temp_degC() << ", Tcore=" << _ht.core_temp_degC() << ", Tbody=" << _ht.body_temp_degC() <<" }" << std::endl;
    _s << "{h_coil_oil="<< htd.h_coil_oil() << ", h_oil_body=" << htd.h_oil_body() << ", h_oil_core=" << htd.h_oil_core() << ", h_body_air=" << htd.h_body_air() << ", h_body_ground=" << htd.h_body_ground() <<" }";
    return _s;
}

std::ostream& Sci::operator << (std::ostream& _s, dt_temp_value const& _t )
{

    _s << "{Tcoil="<< _t.coil_temp_degC() << ", Toil=" << _t.oil_temp_degC() << ", Tcore=" << _t.core_temp_degC() << ", Tbody=" << _t.body_temp_degC() <<" }";
    return _s;
}
/*
std::istream& Sci::operator >> (std::istream& _s, dt_temp_value& _dt )
{
    _s.write ((const char*)&_t, sizeof( _t ) );
    _s.read ((char*)&_dt, sizeof( _dt ) );
    return _s; 
}
*/

void Sci::dt_heat_transfer::rk4( double _dt )
{
    dt_temp_value const& hc_const = hc_;    
    double kcoil1 = _dt * coil_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double koil1  = _dt * oil_ode(  hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double kcore1 = _dt * core_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );
    double kbody1 = _dt * body_ode( hc_const.coil_temp(), hc_const.oil_temp(), hc_const.core_temp(), hc_const.body_temp() );

    double kcoil2 = _dt * coil_ode( hc_const.coil_temp()+kcoil1/2, hc_const.oil_temp()+koil1/2, hc_const.core_temp()+kcore1/2, hc_const.body_temp()+kbody1/2 );
    double koil2  = _dt * oil_ode(  hc_const.coil_temp()+kcoil1/2, hc_const.oil_temp()+koil1/2, hc_const.core_temp()+kcore1/2, hc_const.body_temp()+kbody1/2 );
    double kcore2 = _dt * core_ode( hc_const.coil_temp()+kcoil1/2, hc_const.oil_temp()+koil1/2, hc_const.core_temp()+kcore1/2, hc_const.body_temp()+kbody1/2 );
    double kbody2 = _dt * body_ode( hc_const.coil_temp()+kcoil1/2, hc_const.oil_temp()+koil1/2, hc_const.core_temp()+kcore1/2, hc_const.body_temp()+kbody1/2 );

    double kcoil3 = _dt * coil_ode( hc_const.coil_temp()+kcoil2/2, hc_const.oil_temp()+koil2/2, hc_const.core_temp()+kcore2/2, hc_const.body_temp()+kbody2/2 );
    double koil3  = _dt * oil_ode(  hc_const.coil_temp()+kcoil2/2, hc_const.oil_temp()+koil2/2, hc_const.core_temp()+kcore2/2, hc_const.body_temp()+kbody2/2 );
    double kcore3 = _dt * core_ode( hc_const.coil_temp()+kcoil2/2, hc_const.oil_temp()+koil2/2, hc_const.core_temp()+kcore2/2, hc_const.body_temp()+kbody2/2 );
    double kbody3 = _dt * body_ode( hc_const.coil_temp()+kcoil2/2, hc_const.oil_temp()+koil2/2, hc_const.core_temp()+kcore2/2, hc_const.body_temp()+kbody2/2 );

    double kcoil4 = _dt * coil_ode( hc_const.coil_temp()+kcoil3, hc_const.oil_temp()+koil3, hc_const.core_temp()+kcore3, hc_const.body_temp()+kbody3 );
    double koil4  = _dt * oil_ode(  hc_const.coil_temp()+kcoil3, hc_const.oil_temp()+koil3, hc_const.core_temp()+kcore3, hc_const.body_temp()+kbody3 );
    double kcore4 = _dt * core_ode( hc_const.coil_temp()+kcoil3, hc_const.oil_temp()+koil3, hc_const.core_temp()+kcore3, hc_const.body_temp()+kbody3 );
    double kbody4 = _dt * body_ode( hc_const.coil_temp()+kcoil3, hc_const.oil_temp()+koil3, hc_const.core_temp()+kcore3, hc_const.body_temp()+kbody3 );

    hc_.coil_temp( hc_const.coil_temp() + (kcoil1 + 2*kcoil2 + 2*kcoil3 + kcoil4)/6 );
    hc_.oil_temp ( hc_const.oil_temp()  + (koil1 + 2*koil2 + 2*koil3 + koil4)/6 );
    hc_.core_temp( hc_const.core_temp() + (kcore1 + 2*kcore2 + 2*kcore3 + kcore4)/6 );
    hc_.body_temp( hc_const.body_temp() + (kbody1 + 2*kbody2 + 2*kbody3 + kbody4)/6 );
}



/*
    double kx1 = dt * fx(t, x, y);
    double ky1 = dt * fy(t, x, y);
    
    double kx2 = dt * fx(t+dt/2, x+kx1/2, y+ky1/2 );
    double ky2 = dt * fy(t+dt/2, x+kx1/2, y+ky1/2 );

    double kx3 = dt * fx(t+dt/2, x+kx2/2, y+ky2/2 );
    double ky3 = dt * fy(t+dt/2, x+kx2/2, y+ky2/2 );

    double kx4 = dt * fx(t+dt, x+kx3, y+ky3 );
    double ky4 = dt * fy(t+dt, x+kx3, y+ky3 );


    *x_out = x + (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6;
    *y_out = y + (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6;
*/

	/*
	rk4( dt );
	//calc error
	double x = t0 + dt*i;
	double y2 = std::pow( x * x / 4 + 1, 2);
	std::cout << "STEP: " << x << " " << *this << " error: " << hc_.oil_temp()/y2 - 1 << " y2: " << y2 << std::endl;
	*/
