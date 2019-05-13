#include <iostream>

int main()
{
    for ( int i = 0; i < 5; ++i )
    {
	for (int j = 0; j < 10; ++j )
	{
	    std::cout << "i: " << i << ", j: " << j << std::endl;
	    if ( j == 5 ) break;
	}
    }

    return 0;
}