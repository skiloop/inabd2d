#ifndef _DEFBND_H_
#define _DEFBND_H_
#define  COORDINATES MyStruct
#define MDOUBLE MyDataF
COORDINATES Hzx_xn;
COORDINATES Hzx_xp;

COORDINATES Hzx_yn;
COORDINATES Hzx_yp;

COORDINATES Hzy_xn;
COORDINATES Hzy_xp;

COORDINATES Hzy_yn;
COORDINATES Hzy_yp;


COORDINATES Ceye_xn;
COORDINATES Ceyhz_xn;
COORDINATES Chzxh_xn;
COORDINATES Chzyh_xn;
COORDINATES Chzyex_xn;
COORDINATES Chzxey_xn;

COORDINATES Ceye_xp;
COORDINATES Ceyhz_xp;
COORDINATES Chzxh_xp;
COORDINATES Chzyh_xp;
COORDINATES Chzyex_xp;
COORDINATES Chzxey_xp;

COORDINATES Cexe_yp;
COORDINATES Cexhz_yp;
COORDINATES Chzxh_yp;
COORDINATES Chzyh_yp;
COORDINATES Chzxey_yp;
COORDINATES Chzyex_yp;

COORDINATES Cexe_yn;
COORDINATES Cexhz_yn;
COORDINATES Chzxh_yn;
COORDINATES Chzyh_yn;
COORDINATES Chzxey_yn;
COORDINATES Chzyex_yn;




COORDINATES Ezx_xn;
COORDINATES Ezx_xp;

COORDINATES Ezx_yn;
COORDINATES Ezx_yp;

COORDINATES Ezy_xn;
COORDINATES Ezy_xp;

COORDINATES Ezy_yn;
COORDINATES Ezy_yp;

COORDINATES Cezxe_xn;
COORDINATES Cezxhy_xn;
COORDINATES Chyh_xn;
COORDINATES Chyez_xn;
COORDINATES Cezye_xn;
COORDINATES Cezyhx_xn;

COORDINATES Cezye_yn;
COORDINATES Cezyhx_yn;
COORDINATES Cezxe_yn;
COORDINATES Cezxhy_yn;
COORDINATES Chxh_yn;
COORDINATES Chxez_yn;


COORDINATES Cezxe_xp;
COORDINATES Cezxhy_xp;
COORDINATES Cezye_xp;
COORDINATES Cezyhx_xp;
COORDINATES Chyh_xp;
COORDINATES Chyez_xp;

COORDINATES Cezxe_yp;
COORDINATES Cezxhy_yp;
COORDINATES Cezye_yp;
COORDINATES Cezyhx_yp;
COORDINATES Chxh_yp;
COORDINATES Chxez_yp;

void init_coordinates(COORDINATES *coo, int nx, int ny, MDOUBLE InitValue) {
    int i, j;
    MDOUBLE *p;

    coo->data = (MDOUBLE*) malloc(nx * ny * sizeof (MDOUBLE));
    if (coo->data == NULL) {
        fprintf(stderr, "Alloc space fails!\n");
        exit(0);
    }
    coo->nx = nx;
    coo->ny = ny;
    p = coo->data;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            *p++ = InitValue;
}
#endif

