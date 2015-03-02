#ifndef _STRESS_H
#define _STRESS_H

//	--------------------------------------------------------------------------------------
//	SIGMA11(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis x1
//	SIGMA21(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis x1
//	SIGMA31(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis x1
//	SIGMA12(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis X2
//	SIGMA22(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis X2
//	SIGMA32(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis X2
//	SIGMA13(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis X3
//	SIGMA23(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis X3
//	SIGMA33(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis X3

double  ***sigm11, ***sigm21, ***sigm31,
        ***sigm12, ***sigm22, ***sigm32,
        ***sigm13, ***sigm23, ***sigm33;

#endif 
