#include "gf_cal.h"
#include "bm_decoding.h"
#include "encoding.h"

void main()
{
	unsigned char i = 0;

	systematic_encoding();

	/*transmission through channel*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
	}
	
	bm_decoding();

	return;
}
