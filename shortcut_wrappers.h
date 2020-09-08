#ifndef SHORTCUT_WRAPPERS_H_INCLUDED
#define SHORTCUT_WRAPPERS_H_INCLUDED

//For looping over the system in real space
#define for_xijk \
for(x_i=0;x_i<nx;x_i++) \
	for(x_j=0;x_j<ny;x_j++) \
		for(x_k=0;x_k<nz;x_k++) \
		{ xijk=x_i*nx*ny+x_j*nx+x_k;
#define efor_xijk }
//In the above loop, the order of the looping is changed to give better performance

//For looping over the system in real space
#define fxijk \
for(x_i=0;x_i<nx;x_i++) \
	for(x_j=0;x_j<ny;x_j++) \
		for(x_k=0;x_k<nz;x_k++) \
		{ 
#define efxijk }
//In the above loop, the order of the looping is changed to give better performance

//For looping over the system in k space
#define for_gijk \
for(g_i=0;g_i<nx;g_i++) \
	for(g_j=0;g_j<ny;g_j++) \
		for(g_k=0;g_k<(nz/2+1);g_k++) \
		{ gijk=g_i*ny*(nz/2+1)+g_j*(ny/2+1)+g_k;
#define efor_gijk }

//For looping over the system in k space
#define fgijk \
for(g_i=0;g_i<nx;g_i++) \
	for(g_j=0;g_j<ny;g_j++) \
		for(g_k=0;g_k<(nz/2+1);g_k++) \
		{ 
#define efgijk }

#endif // SHORTCUT_WRAPPERS_H_INCLUDED
