# [3689292520605548544, 723397314333413520, 4912860560424317034, 14720090973493544128, 5413889704333279232, 721988857123348560, 2884037210989275280, 9809174240623242400]
S = []
for X in range(64):
	x=[(X>>i)&1 for i in range(6)]
	g = [0]*8; y = [0]*6

	g[0] = (x[1]) & (x[4]) ^ x[4]
	g[1] = (x[2]^x[5]) & (x[0]^x[1]^x[2]^x[4]^g[0]) 
	g[2] = (x[0]^x[3]^x[4]^x[5]) & (x[0]^x[2]^g[0]) ^ x[2]^ x[4]^ g[1]
	g[3] = (x[1]) & (x[0]^x[3]^x[5]^g[2]) 
	g[4] = (x[4]) & (x[2]^x[5]^g[2]) 
	g[5] = (x[0]^x[2]^x[3]^x[4]) & (x[2]^x[5]) 
	g[6] = (x[0]^x[1]^x[2]^x[3]) & (x[0]^x[5]^g[0]^g[2]^g[5]) 
	g[7] = (x[1]^g[2]) & (x[2]^x[3]^x[5]^g[0]^g[5]) 

	y[0] = g[2] 
	y[1] = x[1] ^ x[4] ^ g[4] ^ g[5] ^ g[6] ^ g[7] 
	y[2] = x[4] ^ x[5] ^ g[0] ^ g[3] ^ g[6] ^ g[7] 
	y[3] = x[0] ^ x[3] ^ x[4] ^ x[5] ^ g[1] ^ g[2] ^ g[3] ^ g[6] 
	y[4] = x[0] ^ x[2] ^ x[4] ^ g[0] ^ g[2] ^ g[3] ^ g[5] ^ g[7] 
	y[5] = g[3] ^ g[4] ^ g[5] ^ g[6] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4)|(y[5]<<5))

print(S)
print([0, 1, 2, 3, 4, 7, 5, 6, 8, 16, 32, 42, 9, 19, 38, 46, 10, 17, 53, 63, 11, 18, 48, 56, 57, 33, 47, 52, 49, 43, 39, 62, 12, 55, 31, 36, 50, 29, 34, 13, 54, 20, 25, 41, 27, 45, 51, 23, 59, 26, 21, 37, 22, 35, 60, 24, 44, 14, 61, 28, 30, 40, 15, 58])
if tuple(S)==tuple([0, 1, 2, 3, 4, 7, 5, 6, 8, 16, 32, 42, 9, 19, 38, 46, 10, 17, 53, 63, 11, 18, 48, 56, 57, 33, 47, 52, 49, 43, 39, 62, 12, 55, 31, 36, 50, 29, 34, 13, 54, 20, 25, 41, 27, 45, 51, 23, 59, 26, 21, 37, 22, 35, 60, 24, 44, 14, 61, 28, 30, 40, 15, 58]):
	print(True)
else:
	print(False)
