#define BOOST_TEST_MODULE individuo
#include <boost/test/included/unit_test.hpp>

#include "individuo.cpp"

BOOST_AUTO_TEST_CASE( individuo_base )
{
	individuo * teste = new individuo(
			157,//x
			255,//y
			7,//especie
			0.1,//taxa_morte
			97,//orientacao
			270,//angvis
			0,//passo = 0, IMPORTANTE para o check
			0.12,//move
			7,//raio
			25,//dens_max
			0.512,//tx_basal
			999// semente
			);
	individuo * copia = new individuo(*teste);
	// testando os construtores
	if ( teste->get_x() != 157 || teste->get_y() != 255) 
		BOOST_ERROR("O construtor do individuo não está posicionando o individuo corretamente");
	if ( copia->get_x() != 157 || copia->get_y() != 255) 
		BOOST_ERROR("O construtor de cópia do individuo não está posicionando o individuo corretamente");

	individuo * naopode;
	BOOST_CHECK_THROW( naopode = new individuo (0,0,5,-1 // taxa de morte negativa
				,0,0,0,0,0,10,0,0), individuo_exception);

	BOOST_CHECK_THROW( naopode = new individuo (0,0,5,0,0,0,0,0,0,0 // densidade maxima zero
				,0,0), individuo_exception);
	
	// testando as propriedades
	BOOST_CHECK( teste->get_tempo() == 0); // tempo ainda nao foi sorteado, deve ser zero
	BOOST_CHECK( copia->get_tempo() == 0);

	teste->set_x(-235.27); teste->set_y(6271.23);
	if (teste->get_x() != -235.27 || teste->get_y() != 6271.23) 
		BOOST_ERROR("Métodos de set_x e set_y não estão posicionando corretamente o indivíduo");

	delete teste, copia;
}

BOOST_AUTO_TEST_CASE( individuo_methods ) 
{
	// testando tempo, todas as taxas = 0
	individuo * teste = new individuo(
			157,255,7,0.0,//taxa_morte
			97,270,4,0.0,//move
			7,25,0.0,//tx_basal
			999);
	teste->update();
	BOOST_CHECK( std::isinf(teste->get_tempo())); // DEVE ter tempo infinito
	delete teste;
	// só taxa basal, habitat BOM
	teste = new individuo(
			157,255,7,0,//taxa_morte
			97,270,4,0,//move
			7,25,0.1,//tx_basal
			999);
	teste->set_habitat(1);
	teste->update();
	BOOST_CHECK(teste->get_tempo() > 0.1 && teste->get_tempo() < 1000);
	delete teste;
	// só taxa basal, habitat RUIM -- deve ser infinito
	teste = new individuo(
			157,255,7,0,//taxa_morte
			97,270,4,0,//move
			7,25,0.1,//tx_basal
			999);
	teste->set_habitat(0);
	teste->update();
	BOOST_CHECK( std::isinf(teste->get_tempo()));
	delete teste;
	// todas as taxas
	teste = new individuo(
			157,255,7,0.1,//taxa_morte
			97,270,4,0.1,//move
			7,25,0.1,//tx_basal
			999);
	teste->update();
	BOOST_CHECK(teste->get_tempo() > 0.1 && teste->get_tempo() < 1000);
}
//EOF
