lp3=function(x,m)
{	
m=paste(m)
switch(m,
'0'=x*0,
'1'=x*0,
'2'=x*0,
'3'=(2*3*20*x^0)*sqrt(7),
'4'=(2*3*4*70*x-2*3*140)*sqrt(9),
'5'=(3*4*5*252*x^2-2*3*4*630*x+2*3*560)*sqrt(11),
'6'=(4*5*6*924*x^3-3*4*5*2772*x^2+2*3*4*3150*x-2*3*1680)*sqrt(13),
'7'=(5*6*7*3432*x^4-4*5*6*12012*x^3+3*4*5*16632*x^2-2*3*4*11550*x+2*3*4200)*sqrt(15),
'8'=(6*7*8*12870*x^5-5*6*7*51480*x^4+4*5*6*84084*x^3-3*4*5*72072*x^2+2*3*4*34650*x-2*3*9240)*sqrt(17),
'9'=(7*8*9*48620*x^6-6*7*8*218790*x^5+5*6*7*411840*x^4-4*5*6*420420*x^3+3*4*5*252252*x^2-2*3*4*90090*x+2*3*18480)*sqrt(19),
'10'=(8*9*10*184756*x^7-7*8*9*923780*x^6+6*7*8*1969110*x^5-5*6*7*2333760*x^4+4*5*6*1681680*x^3-3*4*5*756756*x^2+2*3*4*210210*x-2*3*34320)*sqrt(21)
)
}
