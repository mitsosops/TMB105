import argparse
from pathlib import Path
from collections import Counter
from array import array
import numpy as np
from math import ceil


# the inverse linear interpolation can be used to normalize values of a range to values between 0.0 and 1.0
def inv_lerp(range_min, range_max, value):
    if range_min == range_max:
        return None
    return (value - range_min) / (range_max - range_min)


def compare_to_range(bin_nominator, bin_denominator, low, high, scale):
    if low * bin_denominator <= bin_nominator * scale <= high * bin_denominator:
        return 0
    elif bin_nominator * scale < low * bin_denominator:
        return -1
    elif bin_nominator * scale > high * bin_denominator: 
        return 1


def encode_freqs(freqs: dict):
    byte_array = array('B')
    # add the null terminator to the coded freqs byte array to indicate its end

    max_freq = max(freqs.values())
    num_bytes_for_freq = 0

    while max_freq > 0:
        max_freq >>= 8
        num_bytes_for_freq += 1

    for s, f in freqs.items():
        if s == "\0":
            break

        byte_array.append(ord(s))

        freq_bytes = []
        for _ in range(0, num_bytes_for_freq):
            freq_bytes.append(f & 0xFF)
            f >>= 8

        for b in freq_bytes[::-1]:
            byte_array.append(b)

    byte_array.insert(0, num_bytes_for_freq)

    byte_array.append(0)

    enc_freqs = byte_array.tobytes()

    return enc_freqs


def decode_freqs(enc_freqs: bytes):
    byte_array = list(enc_freqs)
    num_bytes_for_freq = byte_array[0]
    symbol_freqs = {}

    for i in range(1, len(byte_array), num_bytes_for_freq + 1):
        symbol = chr(byte_array[i])
        if symbol == "\0":
            symbol_freqs[symbol] = 1
            break
        freq_barray = []
        for j in range(0, num_bytes_for_freq):
            freq_barray.append(byte_array[i + j + 1])
        freq = int.from_bytes(freq_barray, byteorder='big', signed=False)
        symbol_freqs[symbol] = freq

    return symbol_freqs


def integer_encode(ascii_text: str, symbol_freqs: dict, freq_table: dict, max_freq: int):
    low = 0
    width = 1
    scale = 1

    for s in ascii_text:
        s_width = symbol_freqs[s]
        s_low = freq_table[s] - s_width

        low *= max_freq  # nominator
        scale *= max_freq  # denominator
        low += width * s_low
        width *= s_width  # low + width = high
        # print((s, s_low, s_width))
        # print((low, width, scale))
    return scale, low, low + width


def binary_encode(scale, low, high):
    bin_nominator = 1
    bin_denominator = 1

    # variables that help to keep the analogy of the decimal and binary scale, for very fast comparisons (no multiplications needed as they change along with thir binary counterparts)
    lbd = low # is essentialy low * bin_denominator for each of the following loop
    hbd = high # is essentialy high * bin_denominator for each of the following loop
    sbn = scale # is essentialy scale * bin_nominator for each of the following loop

    cnt_shift = 0
    while True:
        # print("{0:b}".format(bin_nominator))
        # cmp = compare_to_range(bin_nominator, bin_denominator, low, high, scale)
        if lbd <= sbn <= hbd:
            break
        elif sbn < lbd:
            bin_nominator += 1
            bin_nominator <<= 1
            sbn += scale
            sbn <<= 1
            bin_denominator <<= 1
            lbd <<= 1
            hbd <<= 1
            cnt_shift += 1
        elif sbn > hbd:
            bin_nominator -= 1
            sbn -= scale

    print(cnt_shift)
    return bin_nominator, bin_denominator


def encode(ascii_text: str):
    # Create the frequencies dict before adding the null terminator to the text.
    # Append the null terminator to the text and then to the dictionary with a value of 1,
    # to ensure that it is always the last symbol.
    # The null terminator should be always last, to make the frequency table encoding easier
    symbol_freqs = {item[0]: item[1] for item in Counter(ascii_text).most_common()}
    ascii_text += "\0"
    symbol_freqs["\0"] = 1
    # symbol_ids = {item[0]:idx for idx,item in enumerate(symbol_freqs)}
    symbol_count = len(symbol_freqs)
    freq_cumsum = np.cumsum(list(symbol_freqs.values()))
    max_freq = int(max(freq_cumsum))
    freq_table = {list(symbol_freqs.keys())[i]: int(freq_cumsum[i]) for i in range(0, symbol_count)}

    # text_symbols_map = list(map(lambda c: symbol_ids[c],ascii_text))
    # cf = stats.cumfreq(text_symbols_map, numbins= symbol_count)
    # max_freq = int(max(cf.cumcount))
    # freq_table = {item[0]:int(cf.cumcount[item[1]]) for item in symbol_ids.items()}
    
    print(freq_table)
    num_shifts = [ceil(max_freq/(symbol_freqs[s]*2)) for s in symbol_freqs]
    print(num_shifts)

    scale, low, high = integer_encode(ascii_text, symbol_freqs, freq_table, max_freq)

    byte_array = array('B')
    
    # print(bin_nominator)
    # print(bin_denominator)

    # with open('n_d.txt', 'w') as f:
    #     f.write("{0}".format(low))
    #     f.write("\r\n")
    #     f.write("{0}".format(scale))

    # with open('n_d-b.txt', 'w') as f:
    #     f.write("{0:b}".format(low))
    #     f.write("\r\n")
    #     f.write("{0:b}".format(scale))

    bin_nominator, bin_denominator = binary_encode(scale, low, high)

    bin_compact = bin_nominator | bin_denominator
    while bin_compact > 0:
        lsb = bin_compact & 0xFF
        byte_array.insert(0, lsb)
        bin_compact >>= 8

    enc_data = byte_array.tobytes()

    enc_freqs = encode_freqs(symbol_freqs)

    return enc_data, enc_freqs


def decode(enc_data, enc_freqs):
    symbol_freqs = decode_freqs(enc_freqs)
    symbol_count = len(symbol_freqs)
    freq_cumsum = np.cumsum(list(symbol_freqs.values()))
    max_freq = int(max(freq_cumsum))
    freq_table = {list(symbol_freqs.keys())[i]: int(freq_cumsum[i]) for i in range(0, symbol_count)}
    
    print(freq_table)

    bin_nominator = int.from_bytes(enc_data, byteorder='big', signed=False)
    bin_denominator = 1
    while bin_nominator > bin_denominator:
        bin_denominator <<= 1
    bin_denominator >>= 1
    bin_nominator = ~bin_denominator & bin_nominator

    # print(bin_nominator)
    # print(bin_denominator)

    low = 0
    width = 1
    scale = 1
    dec_text = ''
    current_symbol = ''
    while True:
        for s in symbol_freqs:
            s_width = symbol_freqs[s]
            s_low = freq_table[s] - s_width
            
            prev_low = low
            prev_width = width
            prev_scale = scale

            low *= max_freq  # nominator
            low += width * s_low
            width *= s_width  # low + width = high
            scale *= max_freq  # denominator

            if compare_to_range(bin_nominator, bin_denominator, low, low + width, scale) == 0:
                current_symbol = s
                dec_text += s
                break
            else:
                low = prev_low
                width = prev_width
                scale = prev_scale

        if current_symbol == "\0":
            break
    return dec_text


def encode_file(source_file: Path, dest_file: Path = None):
    source_text = ''
    with open(str(source_file), 'r') as f:
        if f.mode == 'r':
            source_text = f.read()

    if source_text == '':
        return None

    enc_data, enc_freqs = encode(source_text)

    if enc_data is None:
        return None

    if dest_file is not None:
        wf = open(str(dest_file), 'wb')
        wf.write(enc_freqs)
        wf.write(enc_data)
        wf.close()

    dec_text = decode(enc_data, enc_freqs)
    
    print(source_text)
    print(dec_text)

    return enc_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Encodes an ascii file using arithmetic coding")
    parser.add_argument('source_file', metavar='<source file>', help='The source file to encode', type=Path)
    parser.add_argument('--out', metavar='<encoded file>',
                        help='The output path (default: <source file>.enc)', type=Path, dest='out_file')

    args = parser.parse_args()
    args.out_file = args.source_file.with_suffix('.enc') if args.out_file is None else args.out_file
    
    encode_file(args.source_file, args.out_file)
