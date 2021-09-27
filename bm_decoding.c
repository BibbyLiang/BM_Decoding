#include <stdio.h>
#include <string.h>
#include "gf_cal.h"
#include "bm_decoding.h"

unsigned char received_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char output_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char syndrome_cal(unsigned char *recv, unsigned char *synd,
								unsigned int cw_len, unsigned int msg_len)
{
	unsigned int i = 0, j = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;

	for(i = 0; i < cw_len - msg_len; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < cw_len; j++)
		{
			tmp = gf_multp(recv[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//printf("%x %x\n", tmp, tmp_sum);
		}
		synd[i] = tmp_sum;
	}
	printf("Syndrome Polynomial:\n");
	for(i = 0; i < cw_len - msg_len; i++)
	{
		printf("%x ", synd[i]);
	}
	printf("\n");
}

int bm_decoding()
{
	unsigned char i = 0, j = 0, tmp = 0xFF, tmp_sum = 0xFF, lambda_root = 0;
	unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char lambda[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char omega[(CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)];
	unsigned char omega_tmp[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1][(CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN - 1)];
	unsigned char err_location[CODEWORD_LEN];
	unsigned char lambda_dev[(CODEWORD_LEN - MESSAGE_LEN) / 2];
	unsigned char err_mag[CODEWORD_LEN];
	unsigned char codeword[CODEWORD_LEN];
	
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(lambda, 0x0, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	memset(omega, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)));
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN); j++)
		{
			omega_tmp[i][j] = 0xFF;
		}
	}
	memset(err_location, 0x0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(lambda_dev, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2));
	memset(err_mag, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(codeword, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);

	unsigned char sigma = 0xFF;
	unsigned char B[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char B_tmp[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char lambda_tmp[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char L = 0;

	memset(B, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	B[1] = 0;
	memset(B_tmp, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	memset(lambda_tmp, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	lambda_tmp[0] = 0;
	memcpy(lambda, lambda_tmp, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));

	printf("Received:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", received_polynomial[i]);
	}
	printf("\n");

	/*compute the syndrome polynomial from received symbols*/
#if 0
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			tmp = gf_multp(received_polynomial[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//printf("%x %x\n", tmp, tmp_sum);
		}
		syndrome[i] = tmp_sum;
	}
	printf("Syndrome Polynomial:\n");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		printf("%x ", syndrome[i]);
	}
	printf("\n");
#else
	syndrome_cal(received_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
#endif
	
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		tmp = 0xFF;
		for(j = 1; j <= L; j++)
		{
			//printf("%x %x %x\n", tmp, syndrome[i - j], lambda[j]);
			tmp = gf_add(tmp, gf_multp(syndrome[i - j], lambda[j]));
		}
		sigma = gf_add(syndrome[i], tmp);
		//printf("%x %x %x\n", sigma, syndrome[i], tmp);

		if(0 != sigma)
		{
			for(j = 0; j < ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1); j++)
			{
				tmp = gf_multp(B[j], sigma);
				lambda[j] = gf_add(lambda[j], tmp);
			}
		}

		if(2 * L < (i + 1))
		{
			L = (i + 1) - L;
			for(j = 0; j < ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1); j++)
			{
				B[j] = gf_div(lambda_tmp[j], sigma);
				//printf("%x %x %x\n", lambda_tmp[j], sigma, B[j]);
			}
		}

		for(j = ((CODEWORD_LEN - MESSAGE_LEN) / 2); j >= 1; j--)
		{
			//printf("%x %x\n", B[j], B[j - 1]);
			B[j] = B[j - 1];
		}
		B[0] = 0xFF;

		memcpy(lambda_tmp, lambda, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
#if 0
		printf("Lambda Polynomial:\n");
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; j++)
		{
			printf("%x ", lambda[j]);
		}
		printf("\n");
		printf("B Polynomial:\n");
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; j++)
		{
			printf("%x ", B[j]);
		}
		printf("\n");
		printf("Sigma: %d\n", sigma);
		printf("L: %d\n", L);
#endif		
	}

	printf("Lambda Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		printf("%x ", lambda[i]);
	}
	printf("\n");

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN; j++)
		{
			omega_tmp[i][j] = gf_multp(lambda[i], syndrome[j]);
			//printf("%d %d: %x %x %x\n", i, j, lambda[i], syndrome[j], omega_tmp[i][j]);
		}
	}
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN; j++)
		{
			//printf("%d %d: %x %x\n", i, j, omega[i + j], omega_tmp[i][j]);
			omega[i + j] = gf_add(omega[i + j], omega_tmp[i][j]);
		}
	}
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		omega[i] = gf_multp(omega[i], gf_mod_single_term(i, CODEWORD_LEN - MESSAGE_LEN));/*there may be some err about remainder calculation*/
	}
	printf("Omega Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		printf("%x ", omega[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		lambda_root = 0xFF;
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; j++)
		{
			//printf("%d %d: %x %x %x %x\n", i, j, lambda_root, lambda[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
			lambda_root = gf_add(lambda_root, gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
		}
		//printf("%d: %x %x\n", i, lambda_root, power_polynomial_table[i + 1][0]);
		if(0xFF == lambda_root)
		{
			err_location[GF_FIELD - 1 - i] = 1;/*inversion*/
		}
	}
	printf("Error Location:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_location[i]);
	}
	printf("\n");

	/*derivative of lambda*/
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		if(0 != ((i + 1) % 2))
		{
			lambda_dev[i] = lambda[i + 1];
		}
	}
	printf("Lambda Derivative Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		printf("%x ", lambda_dev[i]);
	}
	printf("\n");

	/*magnitude of error*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == err_location[i])
		{
			continue;
		}
		else
		{
			tmp = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN - 1) / 2 + 1; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp, omega[j], j * power_polynomial_table[i + 1][0], gf_multp(omega[j], j * power_polynomial_table[i + 1][0]));
				tmp = gf_add(tmp, gf_div(omega[j], j * power_polynomial_table[i + 1][0]));
			}
			tmp_sum = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp_sum, lambda_dev[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
				tmp_sum = gf_add(tmp_sum, gf_div(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
			}
			//printf("%d: %x %x\n", i, tmp, tmp_sum);
			err_mag[i] = gf_div(tmp, tmp_sum);
		}
	}
	printf("Error Magnitude:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_mag[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		codeword[i] = gf_add(codeword[i], err_mag[i]);
	}
	printf("Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", codeword[i]);
	}
	printf("\n");

	return 0;
}
