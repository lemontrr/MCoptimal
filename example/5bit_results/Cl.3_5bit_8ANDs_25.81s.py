# [106956384, 2425751040, 892534784, 4194540838, 2403729792, 2349121836, 209754112, 2181181664]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*8; y = [0]*5

	g[0] = (x[0]^x[1]) & (x[2]^x[3]) 
	g[1] = (x[0]^x[1]^x[2]) & (x[4]^g[0]) 
	g[2] = (x[4]) & (x[0]^x[2]^g[1]) ^ x[4]
	g[3] = (x[0]^g[0]^g[1]) & (x[4]^g[0]^g[2]) ^ x[0]^ x[1]^ x[3]^ x[4]^ g[1]
	g[4] = (x[2]^x[3]^g[1]) & (x[0]^x[4]^g[3]) 
	g[5] = (x[1]^x[2]^x[3]^x[4]^g[1]) & (x[1]^x[4]^g[0]^g[2]^g[4]) 
	g[6] = (x[0]^x[1]^x[3]^g[0]^g[1]) & (x[3]^x[4]^g[0]^g[2]^g[4]) 
	g[7] = (x[3]^g[0]^g[1]^g[2]^g[4]) & (x[1]^x[3]^x[4]^g[1]^g[3]^g[6]) 

	y[0] = x[0] ^ x[1] ^ x[4] ^ g[1] ^ g[2] ^ g[5] ^ g[6] 
	y[1] = x[0] ^ x[1] ^ x[4] ^ g[1] ^ g[2] ^ g[3] ^ g[5] ^ g[6] 
	y[2] = x[2] ^ g[4] ^ g[6] 
	y[3] = x[3] ^ g[1] ^ g[2] ^ g[4] 
	y[4] = x[3] ^ x[4] ^ g[0] ^ g[1] ^ g[2] ^ g[4] ^ g[7] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 6, 7, 8, 5, 16, 18, 29, 9, 15, 28, 26, 10, 30, 20, 19, 23, 31, 24, 11, 12, 22, 27, 17, 13, 21, 14, 25])
if tuple(S)==tuple([0, 1, 2, 3, 4, 6, 7, 8, 5, 16, 18, 29, 9, 15, 28, 26, 10, 30, 20, 19, 23, 31, 24, 11, 12, 22, 27, 17, 13, 21, 14, 25]):
	print(True)
else:
	print(False)
