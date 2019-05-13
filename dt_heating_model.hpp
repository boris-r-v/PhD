#ifndef _Sci_sci_summary_temperature_hpp_
#define _Sci_sci_summary_temperature_hpp_

#include <iostream>
#include <random>

namespace Sci
{
    /**
    *Класс описание конструктивных параметров физических свойств материала дроссель-трансформатора
    * Конструктивные элементы:
    * -coil - обмотка
    * -oil - масло
    * -body - корпус
    * -core - сердечник
    * Характеристика коструктивных элементов
    *	Обмотка: масса(coil_M), теплоемкость(coil_HC), удельное сопротивление(coil_R), площадь поперечного сечения(coil_CSA)[mm^2], длинна(coil_L)[m^2], площадь поверхности контакта с маслом(coil_S)
    *	Сердечник: масса(core_M), теплоемкость(core_HC), площадь поверхности контакта с маслом(core_S)
    *	Масло: масса(oil_M), теплоемкость(oil_HC)
    *	Корпус: масса(body_M), теплоемкость(body_HC), площадь поверхности контакта с маслом(bodу_oil_S), площадь поверхности конвективной теплоотдачи(body_air_S),
    *		площадь поверхности теплоотдачи излучением(body_radiation_S), площадь поверхности кондуктивной теплоотдачи(body_ground_S), площадь нагрева солнцем(body_sun_S)
    *		cтепень серости поверхности корпуса ДТ(body_grayness), номинальный ток ДТ (nominal_curent)	

    */
    class dt
    {
	double coil_M_, coil_HC_, coil_R_, coil_CSA_, coil_L_, coil_S_;
	double core_M_, core_HC_, core_S_;
	double oil_M_, oil_HC_;
	double body_M_, body_HC_, body_oil_S_, body_air_S_, body_radiation_S_, body_ground_S_, body_sun_S_, body_grayness_;
	double nominal_current_;
	std::string title_;
	public:
	    dt()=delete;
	    ~dt()=default;
	    dt( std::string const& title,
		double coil_M, double coil_HC, double coil_R, double coil_CSA, double coil_L, double coil_S,
		double core_M, double core_HC, double core_S,
		double oil_M, double oil_HC,
		double body_M, double body_HC, double body_oil_S, double body_air_S, double body_radiation_S, double body_ground_S, double body_sun_S, double body_grayness,
		double nominal_curent );

	    inline double coil_M() const 	{return coil_M_; }		
	    inline double coil_HC() const 	{return coil_HC_; }		
	    inline double coil_R() const 	{return coil_R_; }		
	    inline double coil_CSA() const 	{return coil_CSA_; }		
	    inline double coil_L() const 	{return coil_L_; }		
	    inline double coil_S() const 	{return coil_S_; }		

	    inline double core_M() const 	{return core_M_; }		
	    inline double core_HC() const 	{return core_HC_; }		
	    inline double core_S() const 	{return core_S_; }		

	    inline double oil_M() const 	{return oil_M_; }		
	    inline double oil_HC() const 	{return oil_HC_; }		

	    inline double body_M() const 	{return body_M_; }		
	    inline double body_HC() const 	{return body_HC_; }		
	    inline double body_oil_S() const 	{return body_oil_S_; }		
	    inline double body_air_S() const 	{return body_air_S_; }		
	    inline double body_radiation_S() const 	{return body_radiation_S_; }		
	    inline double body_ground_S() const 	{return body_ground_S_; }		
	    inline double body_sun_S() const 	{return body_sun_S_; }		
	    inline double body_grayness() const {return body_grayness_; }

	    inline double nominal_current() const {return nominal_current_;}
	    inline std::string title() const {return title_;}
    };
    /**
    * Класс с текущей темпероатурой ДТ
    */
    class dt_temp_value
    {
	    double oil_temp_, coil_temp_, core_temp_, body_temp_;
	public:
	    dt_temp_value();
	    ~dt_temp_value() = default;
	    dt_temp_value( double oil_temp, double coil_temp, double core_temp, double body_temp );

	    inline double oil_temp_degC() const 	{return oil_temp_ - 273; }
	    inline double coil_temp_degC() const 	{return coil_temp_ - 273; }
	    inline double core_temp_degC() const 	{return core_temp_ - 273; }
	    inline double body_temp_degC() const 	{return body_temp_ - 273; }

	    inline double oil_temp() const 	{return oil_temp_; }
	    inline double coil_temp() const 	{return coil_temp_; }
	    inline double core_temp() const 	{return core_temp_; }
	    inline double body_temp() const 	{return body_temp_; }

	    inline void oil_temp(double d) 	{ oil_temp_ = d; }
	    inline void coil_temp(double d) 	{ coil_temp_ = d; }
	    inline void core_temp(double d) 	{ core_temp_ = d; }
	    inline void body_temp(double d) 	{ body_temp_ = d; }

	    /**
	    * Првоверка что температура имеет валидное значение
	    */
	    bool valid() const;

    };

    /**
    * Класс описание параметров, коэффициентов тепловых процессов
    * Конвективная теплоотдача обмотка-масло: 	h_coil_oil
    * Конвективная теплоотдача масло-корпус: 	h_oil_body
    * Конвективная теплоотдача масло-сердечник: h_oil_core
    * Конвективная теплоотдача корпус-воздух:	h_body_air
    * Кондуктивная теплоотдача корпус-грунт:	h_body_ground
    * Поток солнечной радиации:			sun_radiation
    * Температура охлаждающей среды:		external_temp    
    * Облачность, в долях единицы (56% - 0,56): cloud		
    */
    class heat_transfer_data
    {
	    double h_coil_oil_, h_oil_body_, h_oil_core_, h_body_air_, h_body_ground_, sun_radiation_, external_temp_, cloud_, wind_;

	    double get_Pr( double _degC ); // - критерий Прандатля для масла диапазона температур 0-120*С
	    double get_epsilon( double _degK_fluid, double _degK_solid ); //-поправочные коэффициент к теплопредачи

	    std::default_random_engine generator_;
	    std::normal_distribution<double> distribution_;

	public:
	    heat_transfer_data()=delete;
	    ~heat_transfer_data()=default;
	    heat_transfer_data(double h_coil_oil, double h_oil_body, double h_oil_core, double h_body_air, double h_body_ground, double sun_radiation, double external_temp, double cloud, double wind );

	    inline double h_coil_oil() const 	{ return h_coil_oil_; }
	    inline double h_oil_body() const 	{ return h_oil_body_; }
	    inline double h_oil_core() const 	{ return h_oil_core_; }
	    inline double h_body_air() const 	{ return h_body_air_; }
	    inline double h_body_ground() const { return h_body_ground_; }
	    inline double sun_radiation() const { return sun_radiation_; }
	    inline double external_temp() const { return external_temp_; }
	    inline double cloud() const { return cloud_; }
	    inline double wind() const { return wind_; }

	    void h_parametrs_recalc( dt_temp_value const&  );

	    void change_wind( );

    };


    /**
    * Класс реализующий тепловой расчет 
    * начальные условия по температуре часте ДТ: oil_temp_, coil_temp_, core_temp_, bode_temp_
    */
    class dt_heat_transfer
    {
	    dt dt_; 
	    heat_transfer_data htd_;
	    dt_temp_value hc_;
	    double current_;

	    inline double coil_ode	( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp );
	    inline double oil_ode	( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp );
	    inline double core_ode	( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp );
	    inline double body_ode	( double const& _coil_temp, double const& _oil_temp, double const& _core_temp, double const& _body_temp );
	    void rk4 ( double );

	    struct rk45_error{ double coil, oil, core, body; };
	    void rk45 ( double, rk45_error& );
	public:
	    dt_heat_transfer( ) = delete;
	    ~dt_heat_transfer( ) = default;
	    dt_heat_transfer( dt const& dt, heat_transfer_data& htd, dt_temp_value& ic );
	
	    /**
	    *
	    */
	    inline double current() const 	{return current_; }
	    inline void current(double d) 	{ current_ = d; }

	    /**
	    * Метод реализующий расчет
	    * 	current - текущее значение эффективного тока
	    * 	calculate_sec - время рассчета в секундах
    	    */	    
	    void calc( double current, double calculate_sec );
	
	    /**
	    * Получить абсолютные температы элементов
	    */
	    inline double oil_temp() const 		{return hc_.oil_temp(); }
	    inline double coil_temp() const 		{return hc_.coil_temp(); }
	    inline double core_temp() const 		{return hc_.core_temp(); }
	    inline double body_temp() const 		{return hc_.body_temp(); }
	    inline double dt_nominal_current() const 	{return dt_.nominal_current(); }

	    /**
	    * Получить температуру элементов в градусах Цельсия
	    */
	    inline double oil_temp_degC() const 	{return hc_.oil_temp_degC(); }
	    inline double coil_temp_degC() const 	{return hc_.coil_temp_degC(); }
	    inline double core_temp_degC() const 	{return hc_.core_temp_degC(); }
	    inline double body_temp_degC() const 	{return hc_.body_temp_degC(); }
    
	    /**
	    *	Возвращает температуру комопонент ДТ на текущий момент
	    */
	    dt_temp_value get_temp_value() const { return hc_; }
	    /**
	    * Установить температуру компонент
	    */
	    void set_temp_value( dt_temp_value const& _dt );
	    /**
	    * Вернет ссылку на ДТ
	    */
	    dt const& get_dt() const {return dt_; }
	    /**
	    * Вернет ссылку на параметры теплопередачи
	    */
	    heat_transfer_data const& get_heat_transfer_data() const {return htd_; }
	    heat_transfer_data & get_heat_transfer_data() {return htd_; }

    };

    std::ostream& operator << (std::ostream&, dt const& );
    std::ostream& operator << (std::ostream&, heat_transfer_data const& );
    std::ostream& operator << (std::ostream&, dt_heat_transfer const& );
    std::ostream& operator << (std::ostream&, dt_temp_value const& );

}
#endif //_Sci_sci_summary_temperature_hpp_
