#include "sum_temp.hpp"
#include <cmath>
static const double Csb (5.6e-8);

void Sci::heat_transfer_data::h_coil_to_metal_calc( double Toil )
{
    h_coil_oil_ = 0.4019*Toil + 170.09;
    h_oil_body_ = h_coil_oil_;
    h_oil_core_ = h_coil_oil_;   
}

void Sci::heat_transfer_data::h_body_to_air( double Tbody, double Tair )
{
    //h_body_air_=2.1+2.5*std::pow((Tbody - Tair),0.25)*std::pow((Tair/Tbody), 0.25);
    h_body_air_ = 2.1+1.51+1.2834*std::log(Tbody-273);
}


Sci::dt::dt( double coil_M, double coil_HC, double coil_R, double coil_CSA, double coil_L, double coil_S,
		double core_M, double core_HC, double core_S,
		double oil_M, double oil_HC,
		double body_M, double body_HC, double body_oil_S, double body_air_S, double body_radiation_S, double body_ground_S, double body_sun_S, double body_grayness,
		double  nominal_current ):
	coil_M_(coil_M), coil_HC_(coil_HC), coil_R_(coil_R), coil_CSA_(coil_CSA), coil_L_(coil_L), coil_S_(coil_S),
	core_M_(core_M), core_HC_(core_HC), core_S_(core_S),
	oil_M_(oil_M), oil_HC_(oil_HC),
	body_M_(body_M), body_HC_(body_HC), body_oil_S_(body_oil_S), body_air_S_(body_air_S), body_radiation_S_(body_radiation_S), body_ground_S_(body_ground_S), body_sun_S_(body_sun_S), body_grayness_(body_grayness),
	nominal_current_( nominal_current )
{
}


Sci::heat_transfer_data::heat_transfer_data(double h_coil_oil, double h_oil_body, double h_oil_core, double h_body_air, double h_body_ground, double sun_radiation, double external_temp ):
    h_coil_oil_(h_coil_oil),  h_oil_body_(h_oil_body),  h_oil_core_(h_oil_core),  
    h_body_air_(h_body_air),  h_body_ground_(h_body_ground), sun_radiation_(sun_radiation),  
    external_temp_(external_temp)
{
}

Sci::heat_condition::heat_condition( double oil_temp, double coil_temp, double core_temp, double body_temp ):
     oil_temp_(oil_temp),  coil_temp_(coil_temp),  core_temp_(core_temp),  body_temp_(body_temp), current_( 0 )
{
}

Sci::dt_heat_transfer::dt_heat_transfer( dt const& dt, heat_transfer_data& htd, heat_condition& hc ):
    dt_(dt),
    htd_(htd),
    hc_( hc )
{
}

void Sci::dt_heat_transfer::calc( double _current, double _calculating_time )
{
    hc_.current( _current );    
    double dt ( 0.02 );
    double t0 ( 0 );
    for ( double i = t0; i < ( t0 + _calculating_time); i+=dt )
    {
	//std::cout << "STEP: " << i << *this << std::endl;
	rk4( dt );
    }
}

void Sci::dt_heat_transfer::rk4( double _dt )
{
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
    heat_condition const& hc_const = hc_;    
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
    
    htd_.h_coil_to_metal_calc( hc_const.oil_temp() );
    htd_.h_body_to_air( hc_const.body_temp(),  273 );
    //std::cout << " htd: " << htd_ << std::endl;

}

double Sci::dt_heat_transfer::coil_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    //60.0 - считаем время в минутах
    return 60.0*( ( dt_.coil_R()*hc_.current()*hc_.current()*dt_.coil_L()/dt_.coil_CSA() )*(1+0.004*(_coil_temp-293)) - (htd_.h_coil_oil()*dt_.coil_S()*(_coil_temp - _oil_temp) ) ) / ( dt_.coil_M()*dt_.coil_HC() );
}
double Sci::dt_heat_transfer::oil_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    return 60.0*( (htd_.h_coil_oil()*dt_.coil_S()*(_coil_temp - _oil_temp) ) - (htd_.h_oil_body()*dt_.body_oil_S()*(_oil_temp - _body_temp)) - (htd_.h_oil_core()*dt_.core_S()*(_oil_temp - _core_temp)) ) / ( dt_.oil_M()*dt_.oil_HC() );
}
double Sci::dt_heat_transfer::core_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    return 60.0*( htd_.h_oil_core()*dt_.core_S()*(_oil_temp - _core_temp) ) / ( dt_.core_M()*dt_.core_HC() );
}
double Sci::dt_heat_transfer::body_ode( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp )
{
    double colling ( htd_.h_body_air()*dt_.body_air_S()*(_body_temp - htd_.external_temp()) + htd_.h_body_ground()*dt_.body_ground_S()*(_body_temp - htd_.external_temp() ) + 
		dt_.body_grayness()*dt_.body_radiation_S()*Csb*(std::pow(_body_temp,4) - std::pow( htd_.external_temp(),4 ) ) );
    return 60.0*( (htd_.h_oil_body()*dt_.body_oil_S()*(_oil_temp - _body_temp)) + (dt_.body_grayness()*dt_.body_sun_S()*htd_.sun_radiation()) - colling ) / ( dt_.body_M()*dt_.body_HC() );
}


std::ostream& Sci::operator << ( std::ostream& _s, dt const& _d )
{
    return _s;
}
std::ostream& Sci::operator << ( std::ostream& _s, heat_transfer_data const& _h )
{
    _s << "{h_coil_oil: "<< _h.h_coil_oil() << ", h_body_air: " << _h.h_body_air() << ", h_body_ground: " << _h.h_body_ground() << "}";
    return _s;
}

std::ostream& Sci::operator << (std::ostream& _s, dt_heat_transfer const& _ht)
{
//    _s << "{Цельсий: Tbody="<< _ht.body_temp() << ", Tcoil=" << _ht.coil_temp() << ", Tcore=" << _ht.core_temp() << ", Toil=" << _ht.oil_temp() <<" }" << std::endl;
    _s << "{Кельвин: Tbody="<< _ht.body_temp_K() << ", Tcoil=" << _ht.coil_temp_K() << ", Tcore=" << _ht.core_temp_K() << ", Toil=" << _ht.oil_temp_K() <<" }";
//    _s << _ht.body_temp() << "; " << _ht.coil_temp() << "; " << _ht.core_temp() << "; " << _ht.oil_temp() <<"; "<<_ht.body_temp_K() << "; " << _ht.coil_temp_K() << "; " << _ht.core_temp_K() << "; " << _ht.oil_temp_K();
    return _s;
}


