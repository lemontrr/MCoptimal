# [1477787889551294760, 14412117545356309000, 16840579720296941850, 4947804889657918464, 288795545599881924, 648536235702388968, 9908078614235448420, 2702170978839431364, 36714939657758870]
S = []
for X in range(64):
	x=[(X>>i)&1 for i in range(6)]
	g = [0]*9; y = [0]*6

	g[0] = (x[0]^x[1]^x[2]^x[4]) & (x[1]^x[2]^x[3]^x[5]) ^ x[1]^ x[2]^ x[3]^ x[5]
	g[1] = (x[1]^x[3]^x[4]) & (x[1]^x[5]^g[0]) ^ x[1]^ x[3]^ x[4]
	g[2] = (x[2]) & (x[0]^x[1]^x[5]) ^ x[0]^ x[2]^ x[4]^ g[0]^ g[1]
	g[3] = (x[0]^x[3]^x[4]^x[5]) & (x[0]^x[2]^x[3]^g[2]) 
	g[4] = (x[1]^x[3]^x[4]) & (x[1]^x[3]^x[5]^g[2]) 
	g[5] = (x[2]^g[2]) & (x[3]^g[0]^g[4]) 
	g[6] = (x[0]^x[1]^x[5]) & (x[3]^g[0]^g[4]) 
	g[7] = (x[1]^x[5]^g[0]) & (x[1]^x[3]^x[5]^g[2]) 
	g[8] = (x[1]^x[3]^x[4]^g[2]) & (x[0]^x[2]^g[0]^g[4]) 

	y[0] = x[1] ^ x[2] ^ x[3] ^ x[4] ^ g[1] ^ g[2] ^ g[6] 
	y[1] = x[3] ^ x[5] ^ g[0] ^ g[1] ^ g[2] ^ g[3] ^ g[6] ^ g[7] ^ g[8] 
	y[2] = x[0] ^ x[3] ^ x[4] ^ g[1] ^ g[3] ^ g[4] ^ g[6] ^ g[7] ^ g[8] 
	y[3] = x[0] ^ x[1] ^ x[2] ^ x[5] ^ g[0] ^ g[1] ^ g[2] ^ g[4] ^ g[5] ^ g[6] ^ g[7] 
	y[4] = x[1] ^ x[3] ^ x[5] ^ g[0] ^ g[3] ^ g[6] 
	y[5] = x[0] ^ x[1] ^ x[2] ^ g[3] ^ g[4] ^ g[7] ^ g[8] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4)|(y[5]<<5))

print(S)
print([0, 1, 2, 3, 4, 8, 5, 16, 6, 24, 7, 9, 29, 23, 32, 35, 10, 22, 57, 62, 63, 55, 56, 50, 17, 25, 34, 33, 51, 54, 11, 28, 12, 27, 53, 38, 43, 49, 46, 41, 18, 30, 52, 44, 21, 13, 48, 37, 20, 26, 45, 60, 61, 39, 15, 19, 58, 36, 31, 14, 59, 40, 42, 47])
if tuple(S)==tuple([0, 1, 2, 3, 4, 8, 5, 16, 6, 24, 7, 9, 29, 23, 32, 35, 10, 22, 57, 62, 63, 55, 56, 50, 17, 25, 34, 33, 51, 54, 11, 28, 12, 27, 53, 38, 43, 49, 46, 41, 18, 30, 52, 44, 21, 13, 48, 37, 20, 26, 45, 60, 61, 39, 15, 19, 58, 36, 31, 14, 59, 40, 42, 47]):
	print(True)
else:
	print(False)
