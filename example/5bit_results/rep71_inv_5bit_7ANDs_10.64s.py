# [1010565120, 2685012640, 1098357910, 15731280, 847288448, 310444896, 549540480]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*7; y = [0]*5

	g[0] = (x[1]^x[2]) & (x[4]) 
	g[1] = (x[0]) & (x[2]^x[3]^x[4]) 
	g[2] = (x[0]^x[3]^x[4]) & (g[0]^g[1]) ^ x[0]^ x[1]^ x[2]^ x[4]^ g[1]
	g[3] = (x[2]^x[3]) & (x[2]^g[1]) 
	g[4] = (x[4]^g[2]) & (x[2]^x[4]^g[0]^g[3]) 
	g[5] = (x[2]^x[3]^g[0]^g[2]) & (x[3]^x[4]^g[1]^g[3]^g[4]) 
	g[6] = (x[2]^x[4]^g[0]^g[3]) & (x[0]^x[2]^g[1]^g[5]) 

	y[0] = x[0] ^ x[2] ^ x[3] ^ g[0] ^ g[1] ^ g[3] ^ g[5] ^ g[6] 
	y[1] = x[0] ^ x[2] ^ g[1] ^ g[2] ^ g[4] ^ g[5] ^ g[6] 
	y[2] = x[3] ^ g[0] ^ g[3] ^ g[4] 
	y[3] = x[2] ^ g[0] ^ g[3] ^ g[4] ^ g[6] 
	y[4] = x[4] ^ g[4] ^ g[6] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 8, 5, 12, 6, 16, 7, 20, 9, 13, 21, 17, 10, 18, 31, 29, 25, 26, 15, 19, 23, 11, 24, 27, 14, 22, 28, 30])
if tuple(S)==tuple([0, 1, 2, 3, 4, 8, 5, 12, 6, 16, 7, 20, 9, 13, 21, 17, 10, 18, 31, 29, 25, 26, 15, 19, 23, 11, 24, 27, 14, 22, 28, 30]):
	print(True)
else:
	print(False)
