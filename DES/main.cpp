#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// Note
// All of the numbers in the permutation tables are 1-indexed
// like most DES references.

//******************************* Static Data Needed *******************************************

// Applied once at the beginning of the algorithm.
static const unsigned char initial_permutation[64] = {
    58, 50, 42, 34, 26, 18, 10, 2, 60, 52, 44,
    36, 28, 20, 12, 4, 62, 54, 46, 38, 30, 22,
    14, 6, 64, 56, 48, 40, 32, 24, 16, 8, 57,
    49, 41, 33, 25, 17, 9, 1, 59, 51, 43, 35,
    27, 19, 11, 3, 61, 53, 45, 37, 29, 21, 13,
    5, 63, 55, 47, 39, 31, 23, 15, 7};

// Applied once at the end of the algorithm to eliminate effect of initial permutation
static const unsigned char final_permutation[64] = {
    40, 8, 48, 16, 56, 24, 64, 32, 39, 7, 47,
    15, 55, 23, 63, 31, 38, 6, 46, 14, 54, 22,
    62, 30, 37, 5, 45, 13, 53, 21, 61, 29, 36,
    4, 44, 12, 52, 20, 60, 28, 35, 3, 43, 11,
    51, 19, 59, 27, 34, 2, 42, 10, 50, 18, 58,
    26, 33, 1, 41, 9, 49, 17, 57, 25};

// Applied to the half-block at the beginning of the Fiestel function.
static const unsigned char expansion_permutation[48] = {
    32, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9,
    8, 9, 10, 11, 12, 13, 12, 13, 14, 15, 16, 17,
    16, 17, 18, 19, 20, 21, 20, 21, 22, 23, 24, 25,
    24, 25, 26, 27, 28, 29, 28, 29, 30, 31, 32, 1};

// Applied at the end of the Fiestel function.
static const unsigned char fiestel_end_permutation[32] = {
    16, 7, 20, 21, 29, 12, 28, 17, 1, 15, 23,
    26, 5, 18, 31, 10, 2, 8, 24, 14, 32, 27,
    3, 9, 19, 13, 30, 6, 22, 11, 4, 25};

// Converts from full 64-bit key to two key halves: left and right.
// Only 48 bits from the original key are used.
static const unsigned char permuted_choice_1[56] = {
    57, 49, 41, 33, 25, 17, 9, 1, 58, 50, 42, 34,
    26, 18, 10, 2, 59, 51, 43, 35, 27, 19, 11, 3,
    60, 52, 44, 36, 63, 55, 47, 39, 31, 23, 15, 7,
    62, 54, 46, 38, 30, 22, 14, 6, 61, 53, 45, 37,
    29, 21, 13, 5, 28, 20, 12, 4};

// Converts the shifted right and left key halves (concatenated together) into
// the subkey for the round (input into Fiestel function).
static const unsigned char permuted_choice_2[48] = {
    14, 17, 11, 24, 1, 5, 3, 28,
    15, 6, 21, 10, 23, 19, 12, 4,
    26, 8, 16, 7, 27, 20, 13, 2,
    41, 52, 31, 37, 47, 55, 30, 40,
    51, 45, 33, 48, 44, 49, 39, 56,
    34, 53, 46, 42, 50, 36, 29, 32};

// S-Boxes
// Each value represents 4 bits that the 6-bit input is mapped to.
int sbox_1[4][16] =
    {
        14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7,
        0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8,
        4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0,
        15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13};

int sbox_2[4][16] =
    {
        15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10,
        3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5,
        0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15,
        13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9};

int sbox_3[4][16] =
    {
        10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8,
        13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1,
        13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7,
        1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12};

int sbox_4[4][16] =
    {
        7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15,
        13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9,
        10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4,
        3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14};

int sbox_5[4][16] =
    {
        2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9,
        14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6,
        4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14,
        11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3};

int sbox_6[4][16] =
    {
        12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11,
        10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8,
        9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6,
        4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13};

int sbox_7[4][16] =
    {
        4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1,
        13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6,
        1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2,
        6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12};

int sbox_8[4][16] =
    {
        13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7,
        1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2,
        7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8,
        2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11};

// How much the left and right key halves are shifted every round.
static const unsigned char key_shift_amounts[16] = {1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1};

//******************************* End of Static Data *******************************************

/*
 * helper method to make all types of permutations
 *
 * @param input: pointer to array of bytes to apply permutation on it
 * @param table: the table that will applied on input
 * @param output : pointer to array of bytes after permutation
 * @param nbytes : number of bytes of input
 */
void permute(const unsigned char *input, const unsigned char *table, unsigned char *output, unsigned char nbytes)
{
    int idx = 0;
    for (unsigned char i = 0; i < nbytes; i++)
    {
        unsigned char result_byte = 0x00;
        for (unsigned char j = 0; j < 8; j++)
        {
            int val = (*(table + idx)) - 1;
            unsigned char bit_pos_row = val / 8;
            unsigned char bit_pos_col = val % 8;
            unsigned char mask = 0x80 >> bit_pos_col;
            unsigned char result_bit = (input[bit_pos_row] & mask) << bit_pos_col;
            result_byte |= result_bit >> j;
            idx++;
        }
        output[i] = result_byte;
    }
}

/*
 * Make XOR between two byte arrays
 *
 * @param block_a: first byte array
 * @param block_b: second byte array
 * @param output : pointer to array of  bytes after XOR
 * @param nbytes : number of bytes
 */
void XOR(const unsigned char *block_a, const unsigned char *block_b, unsigned char *output, unsigned char nbytes)
{
    for (unsigned char i = 0; i < nbytes; i++)
    {
        output[i] = block_a[i] ^ block_b[i];
    }
}

void des_key_shift(unsigned char key[7], unsigned char output[7], unsigned char amount)
{
    unsigned char mask;

    // Shift bytes circularly.
    // Last byte will be handled specially.
    for (unsigned char i = 0; i < 7; i++)
    {
        output[i] = (key[i] << amount) | (key[i + 1] >> (8 - amount));
    }

    // Middle byte straddles the two key halves.  Set the right 1 or 2 bits of
    // the left key.
    if (amount == 1)
    {
        mask = 0xEF; // 1110 1111
    }
    else
    {
        mask = 0xCF; // 1100 1111
    }
    output[3] &= mask;                             // Zero the bits out
    output[3] |= (key[0] >> (4 - amount)) & ~mask; // Add left 1 or 2 bits

    // Last bit must borrow from left side of right key.
    if (amount == 1)
    {
        mask = 0x01;
    }
    else
    {
        mask = 0x03;
    }
    output[6] = (key[6] << amount) | ((key[3] >> (4 - amount)) & mask);
}

unsigned char get_row(unsigned char a, unsigned char b)
{
    unsigned char result = 0;
    if (a == 0 && b == 0)
    {
        result = 0;
    }
    else if (a == 0 && b == 1)
    {
        result = 1;
    }
    else if (a == 1 && b == 0)
    {
        result = 2;
    }
    else
    {
        result = 3;
    }
    return result;
}
void des_substitution_box(const unsigned char input[6], unsigned char output[4])
{
    unsigned char input_byte, row = 0, col = 0;

    // S-Box 1
    input_byte = (input[0] & 0xFC) >> 2;
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = (input_byte & 0x1E) >> 1;
    output[0] = sbox_1[row][col] << 4;

    // S-Box 2
    input_byte = ((input[0] & 0x03) << 4) | ((input[1] & 0xF0) >> 4);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[0] = output[0] | sbox_2[row][col];

    // S-Box 3
    input_byte = ((input[1] & 0x0F) << 2) | ((input[2] & 0xC0) >> 6);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[1] = sbox_3[row][col] << 4;

    // S-Box 4
    input_byte = (input[2] & 0x3F);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[1] = output[1] | sbox_4[row][col];

    // S-Box 5
    input_byte = (input[3] & 0xFC) >> 2;
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[2] = sbox_5[row][col] << 4;

    // S-Box 6
    input_byte = ((input[3] & 0x03) << 4) | ((input[4] & 0xF0) >> 4);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[2] = output[2] | sbox_6[row][col];

    // S-Box 7
    input_byte = ((input[4] & 0x0F) << 2) | ((input[5] & 0xC0) >> 6);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[3] = sbox_7[row][col] << 4;

    // S-Box 8
    input_byte = (input[5] & 0x3F);
    row = get_row(((input_byte & 0x20) >> 5), (input_byte & 0x01));
    col = ((input_byte & 0x1E) >> 1);
    output[3] = output[3] | sbox_8[row][col];
}

void des_fiestel(const unsigned char input[4], const unsigned char subkey[6], unsigned char output[4])
{
    unsigned char expanded[6];
    unsigned char sbox_output[4];

    // TODO: Can expansion be done faster than a dumb permutation algorithm?
    permute(input, expansion_permutation, expanded, 6);
    // TODO: Can XOR and sbox be combined?
    XOR(expanded, subkey, expanded, 6);
    des_substitution_box(expanded, sbox_output);
    permute(sbox_output, fiestel_end_permutation, output, 4);
}

void generate_keys(unsigned char key[8], unsigned char keys[16][6])
{
    unsigned char key_halves_a[7]; // left key + right key
    unsigned char key_halves_b[7]; // Also left key + right key
    unsigned char subkey[6];
    unsigned char fiestel_output[4];

    // TODO: Pre shift permuted_choice_1 to eliminate left shift to generate
    //       the first subkey.  Or maybe even have a lookup table for each
    //       subkey.
    permute(key, permuted_choice_1, key_halves_a, 7);
    for (unsigned char i = 0; i < 16; i += 2)
    {
        // Generate key (even round)
        des_key_shift(key_halves_a, key_halves_b, key_shift_amounts[i]);
        permute(key_halves_b, permuted_choice_2, subkey, 6);
        for (size_t j = 0; j < 6; j++)
        {
            keys[i][j] = subkey[j];
        }

        // Generate key (odd round)
        des_key_shift(key_halves_b, key_halves_a, key_shift_amounts[i + 1]);
        permute(key_halves_a, permuted_choice_2, subkey, 6);
        for (size_t j = 0; j < 6; j++)
        {
            keys[i + 1][j] = subkey[j];
        }
    }
}

void des_encrypt(unsigned char block[8], unsigned char keys[16][6], unsigned char output[8])
{
    // TODO: This whole program could probably benifit from using larger
    //       datatypes that char.
    unsigned char fiestel_output[4];

    // left_block and right_block must be beside each other in memory, so the
    // memory is allocated in one chunk and the left and right block use 4
    // bytes each.
    unsigned char *left_block = new unsigned char[8];
    unsigned char *right_block = &left_block[4];
    unsigned char blockOut[8];

    permute(block, initial_permutation, blockOut, 8);
    for (size_t i = 0; i < 4; i++)
    {
        left_block[i] = blockOut[i];
        right_block[i] = blockOut[i + 4];
    }

    // Calculate 16 Rounds
    // Each loop iteration calculates two rounds.  This way there are no
    // memcoppies at the end of each round to for example switch right_block
    // and left_block.
    for (unsigned char i = 0; i < 16; i += 2)
    {

        // Round calculation (even round)
        des_fiestel(right_block, keys[i], fiestel_output);
        XOR(fiestel_output, left_block, left_block, 4);

        // Round calculation (odd round)
        des_fiestel(left_block, keys[i + 1], fiestel_output);
        XOR(fiestel_output, right_block, right_block, 4);
    }

    // swap left and right
    for (size_t i = 0; i < 4; i++)
    {
        unsigned char temp = left_block[i];
        left_block[i] = left_block[i + 4];
        left_block[i + 4] = temp;
    }

    permute(left_block, final_permutation, output, 8);
}

void des_decrypt(unsigned char block[8], unsigned char keys[16][6], unsigned char output[8])
{
    // TODO: This whole program could probably benifit from using larger
    //       datatypes that char.
    unsigned char fiestel_output[4];

    // left_block and right_block must be beside each other in memory, so the
    // memory is allocated in one chunk and the left and right block use 4
    // bytes each.
    unsigned char *left_block = new unsigned char[8];
    unsigned char *right_block = &left_block[4];
    unsigned char blockOut[8];

    permute(block, initial_permutation, blockOut, 8);

    for (size_t i = 0; i < 4; i++)
    {
        left_block[i] = blockOut[i];
        right_block[i] = blockOut[i + 4];
    }

    // Calculate 16 Rounds
    // Each loop iteration calculates two rounds.  This way there are no
    // memcoppies at the end of each round to for example switch right_block
    // and left_block.
    for (unsigned char i = 0; i < 16; i += 2)
    {

        // Round calculation (even round)
        des_fiestel(right_block, keys[16 - i - 1], fiestel_output);
        XOR(fiestel_output, left_block, left_block, 4);

        // Round calculation (odd round)
        des_fiestel(left_block, keys[16 - i - 2], fiestel_output);
        XOR(fiestel_output, right_block, right_block, 4);
    }

    // swap left and right
    for (size_t i = 0; i < 4; i++)
    {
        unsigned char temp = left_block[i];
        left_block[i] = left_block[i + 4];
        left_block[i + 4] = temp;
    }

    permute(left_block, final_permutation, output, 8);
}

void read_bytes_hex(FILE *file, unsigned char *text, int nbytes)
{
    for (size_t i = 0; i < nbytes; i++)
    {
        fscanf(file, "%02x", &text[i]);
    }
}

bool read_bytes_bin(FILE *file, unsigned char *text, int nbytes)
{
    for (size_t i = 0; i < nbytes; i++)
    {
        if (fscanf(file, "%c", &text[i]) == EOF)
        {
            return 0;
        }
    }
    return 1;
}

int main(int argc, char *argv[])
{

    // true for encryption / false for decryption
    bool mode = true;

    if (argv[1] == NULL || argv[2] == NULL || argv[3] == NULL || argv[4] == NULL)
    {
        printf("Error: input are not correct");
        printf("\nExample usage: DES.exe (encrypt|decrypt) [TEXT_FILE KEY_FILE OUT_FILE]\n\n");
        return -1;
    }

    if (strcmp(argv[1], "encrypt") == 0)
    {
        mode = true;
    }
    else if (strcmp(argv[1], "decrypt") == 0)
    {
        mode = false;
    }
    else
    {
        printf("Error: input is not correct");
        printf("\nExample usage: DES (encrypt|decrypt) [TEXT_FILE KEY_FILE OUT_FILE]\n\n");
        return -1;
    }

    FILE *input_key = fopen(argv[3], "rb");
    FILE *input_text = fopen(argv[2], "rb");
    FILE *output_text = fopen(argv[4], "wb");

    unsigned char plaintext[8];
    unsigned char key[8];
    unsigned char keys[16][6];
    unsigned char ciphertext[8];

    // fread(key, sizeof key, 1, input_key);
    read_bytes_hex(input_key, key, 8);
    free(input_key);
    generate_keys(key, keys);
    while (read_bytes_bin(input_text, plaintext, 8))
    {
        if (mode)
        {
            des_encrypt(plaintext, keys, ciphertext);
               printf("%s\n",ciphertext);
        }
        else
            des_decrypt(plaintext, keys, ciphertext);

        printf("%s\n",ciphertext);
        fwrite(ciphertext, sizeof ciphertext, 1, output_text);
    }

    fclose(input_key);
    fclose(input_text);
    fclose(output_text);
}
