import math

def calculate_tm(seq, na_conc):
    """
    Calculate the melting temperature (Tm) of a DNA primer using the salt concentration method.
    
    :param seq: str, DNA sequence (e.g., "ATCGATCG")
    :param na_conc: float, sodium ion concentration in mM
    :return: float, melting temperature in Celsius
    """
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')

    tm = 2 * (a_count + t_count) + 4 * (g_count + c_count) + \
         math.log10(na_conc) - 0.3 * na_conc + 16.6

    return tm

# 使用示例
dna_sequence = "AGAGTTTGATCCTGGCTCAG"
sodium_concentration = 50  # 单位为 mM
tm_value = calculate_tm(dna_sequence, sodium_concentration)
print(f"The Tm value for the primer {dna_sequence} is {tm_value:.2f} °C.")