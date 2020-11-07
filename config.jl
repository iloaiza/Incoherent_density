# Unit conversion constants
au2ev = 27.2113845
ev2au = 1/au2ev
au2cm = 219474.63
cm2au = 1/au2cm
ev2cm = 8065.540106923572

# Parameters for Rhodopsin model
const E0rhod = 0.0
const E1rhod = 2.48*ev2au
const V0rhod = 3.6*ev2au
const V1rhod = 1.09*ev2au
const κrhod = 0.1*ev2au
const λrhod = 0.19*ev2au
const Iinvrhod = 4.84e-4*ev2au
const ω2rhod = 0.19*ev2au

# Parameters for 24-D pyrazine model
# Order of elements:
# [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24]
# [6a	1	9a	8a	2	10a	4	5	6b	3	8b	7b	16a	17a	12	18a	19a	13	18b	14	19b	20b	16b	11]
# AiPIR = []

# Quick multiDim parameters for LVC model
N2 = [10,10]; Ai2 = [0,0]; Bi2 = [0,1]; Ω2 = [2,4]; Ci2 = [0,0]; Δ2 = 0.1
N3 = [7,7,7]; Ai3 = [0,0,0]; Bi3 = [0.7,1,1.2]; Ω3 = [2,4,3.1]; Ci3 = [0,0,0]; Δ3 = 0.1
N5 = [8,4,8,4,3] .- 1; Ω5 = [1,1.7,2,2,3]; Ai5 = zeros(5); Bi5 = [3,1,0.5,1,0]; Ci5 = zeros(5); Δ5 = 0.1


DATAFOLDER = "./"