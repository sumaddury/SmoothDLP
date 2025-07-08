#pragma once
#include <gmpxx.h>
#include <cstdint>

mpz_class randSmooth(uint32_t  y, const mpz_class& B);

int legendreSymbol(const mpz_class& a_mp, const mpz_class& p_mp);
