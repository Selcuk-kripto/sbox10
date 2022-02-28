// This code finds difference distribution table (DDT) and linear approximation table (LAT) 
// of some rotation symmetric S-boxes and their concatenations (in dimension 10), which are 
// obtained by the steepest descent like iterative search algorithm.

// Sel�uk Kavut
// Dept. of Computer Engineering
// Bal�kesir University, Turkey

// The code is compiled by Microsoft Visual C++ 2010.

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"

//n: Number of variables
#define n 10
#define N 1024

int *DDT   =(int *)malloc(N*N*sizeof(int));
int *LAT=(int *)malloc(N*N*sizeof(int));
int *T=(int *)malloc(N*sizeof(int));
int *S=(int *)malloc(N*sizeof(int));
int *FW=(int *)malloc(N*sizeof(int));
int *DB=(int *)malloc(N*n*sizeof(int));
int *Sb=(int *)malloc(N*n*sizeof(int));

//Result are written to the files "DDT.txt" and "LAT.txt".
FILE *outddt=fopen("DDT.txt", "w");
FILE *outlat=fopen("LAT.txt", "w");


//S-box (S_1) (rotation symmetric) with nonlinearity 456, absolute indicator 192, algebraic degree 9, and differential uniformity 10 
int RSSB1[]={1023,552,81,993,162,489,963,929,324,456,978,926,903,879,835,802,648,652,912,554,933,288,829,870,783,239,735,351,647,538,581,616,273,660,281,82,801,141,85,241,843,962,576,886,635,370,717,65,543,338,478,977,447,679,702,111,271,748,53,412,139,527,209,512,546,808,297,149,562,257,164,832,579,942,282,764,170,379,482,713,663,20,901,662,129,945,749,45,247,721,740,184,411,225,130,590,63,578,676,825,956,240,931,84,894,892,335,300,381,704,222,794,542,13,473,668,106,1007,824,470,278,911,31,904,418,920,1,968,69,57,593,925,594,376,298,605,101,885,514,218,328,623,641,1009,135,750,861,766,564,468,505,989,340,1013,758,171,964,249,403,472,303,420,40,519,779,495,301,28,258,675,867,744,475,661,90,884,494,971,419,757,457,131,368,559,822,751,450,305,260,820,157,29,126,674,133,125,329,425,627,809,889,943,480,367,839,55,168,811,765,949,761,880,670,201,600,237,762,345,385,603,444,584,565,194,61,551,26,642,946,896,313,693,212,969,991,265,625,972,940,364,556,155,799,700,62,304,785,525,836,707,817,86,2,476,913,818,138,378,114,987,163,72,827,646,165,291,752,604,596,937,187,899,202,320,747,862,5,492,436,312,656,60,223,176,259,1019,995,230,270,94,477,923,699,117,509,318,105,891,936,421,1010,800,955,205,680,362,1003,781,493,306,342,146,905,224,498,243,806,76,944,119,606,18,840,490,80,123,15,1022,535,285,990,200,602,588,56,19,516,798,327,147,711,682,465,922,950,837,299,638,180,619,745,437,988,731,919,393,838,279,491,296,914,173,262,561,736,506,95,927,621,720,479,723,900,430,610,784,520,934,617,678,314,725,58,985,252,883,325,939,266,878,250,563,658,607,850,23,231,549,595,113,755,283,863,405,960,93,734,166,655,941,110,459,336,161,599,778,507,587,875,22,499,692,737,414,317,25,402,847,177,198,474,852,501,898,690,441,770,1001,183,469,888,380,145,683,107,639,388,537,122,438,79,966,52,834,261,413,869,893,769,994,626,846,363,864,424,531,915,137,959,697,530,49,227,394,921,349,857,67,728,369,89,34,310,856,575,667,377,395,124,422,608,1020,547,666,27,66,649,233,391,35,611,854,172,789,4,232,952,204,803,871,613,192,276,1008,756,976,228,463,951,401,326,277,144,435,631,687,269,308,330,41,582,632,481,443,185,544,169,1000,851,567,374,206,775,256,404,586,640,416,471,382,701,868,10,331,984,534,872,92,624,295,289,924,120,42,446,150,352,397,518,334,1015,235,967,452,460,484,540,974,188,814,954,109,823,1016,375,383,234,1006,1018,597,636,236,210,771,759,14,849,372,842,442,997,890,577,791,887,664,410,526,337,574,724,916,983,695,539,917,986,440,612,630,684,813,292,97,787,321,448,858,996,644,486,182,589,350,152,774,865,43,238,409,189,1005,36,323,657,302,980,961,160,431,246,156,30,88,1021,115,47,973,570,159,957,722,400,614,181,902,153,73,112,633,38,571,9,245,573,511,654,100,294,521,399,585,341,461,930,319,821,730,877,708,651,148,598,792,253,975,360,873,215,392,467,339,874,1004,953,981,439,793,815,523,786,568,653,714,558,83,982,741,592,389,805,11,346,207,524,935,99,426,449,732,1012,746,190,853,831,780,219,483,417,718,958,497,423,432,777,580,860,536,197,686,545,696,17,428,845,709,211,510,333,33,628,529,427,906,116,102,947,96,504,488,743,712,650,729,855,154,532,316,733,272,500,795,103,128,293,208,191,434,677,267,46,659,462,21,75,710,167,629,226,242,487,407,566,508,703,503,810,118,897,7,186,221,445,907,332,263,287,458,859,970,220,315,918,560,672,429,322,91,175,387,533,716,1014,673,151,992,727,78,44,569,998,591,361,307,451,548,828,797,634,767,50,772,804,742,671,365,354,74,396,999,948,196,681,502,1002,908,773,284,357,553,882,706,517,615,979,213,366,373,938,390,753,359,760,216,290,268,343,348,214,866,255,528,776,453,51,48,244,356,876,77,158,136,909,64,104,217,645,841,522,355,826,121,715,254,763,59,515,622,965,643,229,485,669,280,726,557,705,358,848,496,39,796,807,665,274,910,895,386,371,694,37,1011,98,251,454,142,788,353,819,618,698,195,691,108,134,174,433,264,738,24,178,550,68,32,620,932,689,572,127,541,311,833,754,140,790,179,248,398,844,455,193,347,1017,637,71,688,309,609,54,87,132,12,275,16,466,286,782,928,70,601,199,739,685,830,344,816,555,6,8,143,464,812,881,415,408,3,583,406,719,513,203,768,384,0};
//S-box (S_2) (rotation symmetric) with nonlinearity 454, absolute indicator 184, algebraic degree 9, and differential uniformity 8 
int RSSB2[]={0,401,802,774,581,572,525,783,139,279,121,152,27,636,543,301,278,478,558,365,242,136,304,287,54,571,249,938,63,503,602,909,556,528,956,423,93,369,730,953,484,2,272,664,608,810,574,1018,108,904,119,149,498,871,853,999,126,180,1006,227,181,848,795,717,89,925,33,568,889,537,846,274,186,389,738,205,437,507,883,158,968,408,4,722,544,970,305,991,193,790,597,49,125,932,1013,733,216,237,785,924,238,620,298,448,996,190,719,289,683,875,975,835,252,292,360,229,989,570,454,261,362,773,673,300,567,1012,411,899,178,930,827,455,66,256,113,534,755,688,51,866,669,791,548,736,372,86,778,426,453,569,410,550,874,339,1014,243,743,546,316,160,913,555,816,403,8,561,421,118,65,460,917,928,610,281,959,794,386,211,557,897,171,44,98,25,250,308,841,843,1003,699,443,940,432,854,474,329,547,173,825,334,476,666,217,120,596,604,896,823,969,951,380,457,415,613,578,436,343,921,727,918,927,757,647,110,504,829,584,14,720,707,458,726,955,276,117,919,908,624,522,554,724,482,523,892,323,430,600,74,111,257,1001,420,822,884,775,402,356,143,837,159,631,34,910,893,132,348,512,714,226,985,45,212,487,390,353,894,102,754,709,233,315,155,559,986,73,654,449,253,744,64,172,965,533,398,852,648,906,396,115,326,820,11,77,942,725,299,678,151,1005,409,486,445,463,944,69,156,632,619,320,221,803,520,87,502,609,700,806,675,16,611,99,770,842,358,236,922,130,175,920,601,811,341,833,9,197,263,562,123,895,36,565,94,772,788,422,41,91,29,771,764,342,202,88,957,196,492,50,270,500,511,616,838,659,144,663,351,983,214,375,799,886,376,857,679,864,19,685,373,948,83,658,412,71,665,346,134,627,164,668,549,952,364,309,414,434,116,240,489,169,15,185,882,769,1010,623,564,915,345,879,379,760,808,914,531,830,352,203,399,133,759,872,995,686,784,819,366,431,945,813,407,831,200,491,187,271,57,220,560,1008,977,635,510,145,1021,28,336,417,418,391,577,916,283,429,644,887,590,552,859,234,576,815,981,793,606,225,473,21,381,85,12,425,863,964,122,23,856,761,337,646,477,860,456,177,127,148,297,222,475,514,737,979,481,840,192,621,359,745,388,527,670,804,518,712,387,286,903,651,76,318,662,239,694,68,655,797,469,763,966,264,723,696,988,1,332,405,509,452,586,947,1011,90,625,424,870,974,284,780,137,706,614,765,79,204,361,485,1007,395,536,466,878,630,462,310,224,95,656,949,929,146,626,285,642,898,150,506,961,465,739,128,267,344,433,907,368,43,213,796,275,681,633,273,80,789,713,792,59,230,464,652,397,617,960,22,524,154,933,861,470,427,676,598,167,333,60,302,923,987,740,818,218,972,459,890,55,926,7,865,363,138,971,312,277,241,446,215,37,640,210,442,201,583,591,17,958,174,357,1004,106,195,447,377,628,589,493,327,638,32,994,199,324,198,163,517,471,661,587,716,734,472,78,821,622,260,251,350,849,817,385,179,461,599,812,682,516,643,573,18,47,394,532,526,382,101,990,246,135,767,419,72,687,107,911,188,851,521,698,553,206,844,67,82,786,182,207,58,756,519,441,505,282,684,701,404,777,176,711,891,1009,392,183,984,715,100,605,540,280,1000,255,1022,168,209,800,653,322,295,941,288,1002,303,748,702,6,943,61,428,680,750,228,575,660,749,880,752,96,691,194,335,259,705,963,38,331,347,839,746,483,873,494,166,766,293,1017,824,435,142,580,307,551,692,1015,268,439,231,112,328,976,313,321,75,992,881,645,728,184,618,649,828,40,868,541,232,710,480,262,978,235,338,595,30,973,370,109,741,539,515,693,997,650,223,530,105,612,807,479,690,53,735,314,758,319,497,162,593,747,805,367,39,311,637,936,704,742,406,258,798,535,266,191,495,579,721,855,967,937,349,103,545,393,615,378,732,141,862,900,867,1016,603,869,814,140,639,84,400,161,982,501,374,3,542,340,114,330,440,48,97,641,993,677,931,753,247,383,1020,729,290,787,1019,731,56,488,672,496,834,92,836,20,782,355,131,629,809,998,566,781,858,325,265,306,751,538,157,671,81,885,695,667,468,371,129,779,607,801,939,980,563,708,189,582,450,508,946,70,42,592,762,513,170,165,24,832,850,888,703,876,905,877,244,248,46,10,689,826,499,902,674,153,269,847,954,845,697,901,912,490,354,291,254,35,296,768,594,416,444,438,950,124,5,413,451,588,935,934,962,245,657,529,384,208,219,62,718,294,467,634,776,104,31,147,317,52,585,26,13,1023};

//S-box (S_3) (concatenation of rotation symmetric S-boxes) with nonlinearity 456, absolute indicator 192, algebraic degree 9, and differential uniformity 12 
int CONC1[]={0,170,340,572,169,263,632,638,338,376,15,176,752,458,764,434,165,201,241,977,30,910,352,294,992,361,405,938,1016,1021,357,968,330,300,402,742,482,187,931,334,60,457,797,975,193,192,77,54,961,615,211,516,299,212,853,410,1009,765,1019,483,203,1001,913,929,149,47,89,109,293,121,972,863,453,292,374,448,839,905,157,525,120,437,403,788,571,736,927,623,386,7,384,90,154,289,108,942,899,573,718,640,422,657,520,829,87,590,424,236,683,817,309,431,995,234,1018,985,1015,95,455,1002,406,616,979,754,803,722,835,412,298,449,94,370,178,867,218,767,75,430,242,48,921,53,703,762,395,158,73,738,237,568,385,200,655,676,787,716,314,407,538,820,240,984,363,141,295,526,553,485,630,771,960,61,831,244,734,600,261,766,14,199,257,465,180,111,308,889,67,490,216,864,861,717,775,22,634,693,924,1017,768,252,333,56,802,139,528,285,635,670,174,6,668,793,337,327,472,189,854,209,611,658,107,444,351,247,967,210,468,345,1013,998,947,862,1007,268,190,486,399,427,981,987,301,353,720,253,935,898,996,958,583,886,932,443,647,821,313,680,85,542,387,575,188,88,229,217,356,1000,711,147,436,725,1022,740,150,627,349,167,484,999,96,27,819,514,106,205,894,497,1012,976,279,310,316,943,146,224,964,774,474,650,624,823,259,45,400,727,798,576,840,926,551,118,920,471,117,1004,303,757,564,633,617,206,480,185,945,895,215,24,282,637,79,625,540,100,594,614,459,666,748,326,519,498,897,286,122,556,639,355,488,311,956,245,688,870,11,858,1020,126,28,325,398,591,3,908,419,350,360,585,222,379,105,428,755,687,134,243,469,1005,432,382,705,735,699,477,922,596,527,799,44,364,756,329,874,626,825,339,1011,269,513,358,504,744,155,983,112,643,581,923,278,875,544,719,59,491,758,890,828,103,348,959,12,830,824,50,563,589,163,249,143,534,433,411,378,691,685,63,418,807,710,175,804,445,214,855,377,1014,191,879,494,554,911,182,420,569,425,390,179,628,1003,833,973,949,871,501,701,307,991,671,25,806,380,523,461,857,287,915,343,478,939,763,951,533,91,796,195,570,928,986,506,409,847,659,773,940,969,239,893,778,654,541,749,460,841,726,375,645,782,230,619,834,115,673,848,511,1023,906,789,559,555,291,606,99,598,562,71,34,700,462,198,560,684,43,612,64,142,275,68,795,888,347,413,296,396,729,608,464,856,354,86,733,712,481,128,595,284,97,39,70,136,707,567,251,753,982,183,318,315,104,81,548,281,517,946,660,704,84,417,184,689,646,197,235,172,76,954,832,912,877,451,1006,256,127,678,332,57,31,194,565,78,826,140,741,272,887,902,721,622,698,502,181,994,4,941,467,366,675,125,792,119,505,208,588,162,812,584,440,51,818,522,904,869,842,808,962,896,101,168,761,323,290,368,723,866,456,780,371,394,452,470,694,344,248,152,944,885,26,641,21,801,19,731,415,391,846,989,814,1,936,254,587,844,850,153,328,114,113,62,262,388,851,618,852,156,980,629,233,280,421,970,783,33,98,751,706,781,151,930,129,732,811,884,65,493,574,362,827,965,132,8,37,859,392,423,786,221,957,838,602,250,777,561,607,238,566,499,713,416,93,664,72,324,651,601,219,656,5,369,507,102,881,613,642,532,686,785,604,715,978,661,40,593,260,901,784,769,439,202,397,336,760,1010,776,135,426,69,479,225,750,934,473,709,791,401,305,537,17,231,536,277,32,393,909,429,148,876,232,177,878,496,809,304,35,865,381,747,159,52,530,770,586,42,92,579,373,38,672,950,759,319,166,271,794,669,882,955,872,605,346,2,489,849,652,508,550,662,220,665,708,677,737,306,892,145,873,228,441,226,603,124,728,13,266,265,463,679,663,724,805,681,164,312,131,937,682,746,372,466,903,49,609,331,320,917,288,543,925,66,274,196,649,990,557,900,815,539,868,302,36,837,365,258,509,952,577,599,558,745,20,130,648,475,454,636,644,213,495,631,492,907,408,264,524,16,966,74,116,695,916,273,446,335,521,549,46,442,592,891,83,653,697,692,173,500,582,531,110,610,880,702,948,476,813,620,133,487,843,914,82,321,597,186,963,816,160,144,974,137,836,790,919,690,18,438,510,800,535,10,580,227,578,503,246,204,518,739,58,714,223,772,23,552,297,860,342,547,55,696,730,918,322,933,41,810,993,80,743,674,971,9,255,779,546,545,123,515,29,367,267,404,171,283,621,161,276,1008,883,997,383,529,317,270,389,341,822,138,953,447,414,450,667,988,207,845,359,435,512};
//S-box (S_4) (concatenation of rotation symmetric S-boxes) with nonlinearity 454, absolute indicator 184, algebraic degree 9, and differential uniformity 8 
int CONC2[]={511,127,254,28,508,661,56,888,505,670,810,799,112,435,753,515,499,544,828,331,597,181,575,16,224,447,359,686,994,862,518,239,487,514,576,332,633,188,151,327,682,672,362,262,638,935,32,780,448,202,383,816,207,33,860,729,965,843,701,886,524,480,478,59,463,915,516,503,640,532,153,745,754,365,376,467,302,267,143,529,852,453,832,349,213,99,13,651,764,159,847,821,64,703,537,369,385,372,404,531,255,437,609,743,414,92,66,938,697,905,946,504,907,122,663,631,890,35,749,626,536,648,449,999,445,684,118,947,415,677,807,492,520,173,495,855,768,47,552,1001,306,136,978,120,996,517,219,450,241,408,423,943,93,237,23,738,286,392,546,555,681,171,395,34,641,102,187,98,426,281,198,18,26,72,790,782,1016,104,318,685,671,288,619,825,128,602,894,441,562,570,227,1021,259,995,233,917,297,416,550,942,510,250,363,870,706,693,974,1020,317,637,184,668,132,130,853,119,882,940,787,950,869,742,497,303,791,1,244,813,814,872,750,583,757,1018,70,819,986,231,740,595,560,712,784,399,387,744,975,818,379,397,856,380,236,1015,871,217,319,14,842,700,591,911,473,769,528,421,346,8,479,599,687,375,513,166,94,419,592,131,979,646,101,664,272,876,933,699,240,285,969,507,522,884,438,489,389,776,482,430,305,837,335,922,863,440,186,777,474,883,46,725,964,252,61,827,273,569,580,1011,598,985,850,246,342,939,279,1012,68,60,770,225,204,983,374,625,196,789,341,17,51,49,396,9,36,647,52,854,144,924,557,476,541,1022,1009,970,208,727,125,691,858,766,830,590,65,315,726,731,627,407,256,918,692,803,765,921,371,809,612,455,628,665,454,190,1019,364,7,606,967,896,466,4,811,443,83,465,321,579,588,694,861,398,509,698,500,644,215,930,717,220,900,953,874,126,925,796,1017,1004,123,981,762,30,368,1003,824,906,264,280,260,835,683,718,238,767,741,875,857,639,551,413,877,459,715,913,972,916,483,844,95,182,559,704,2,477,488,801,603,199,605,578,721,110,988,63,654,758,1002,15,1013,709,140,929,615,895,949,831,462,485,968,714,678,91,608,494,912,355,545,55,287,635,263,866,976,959,927,498,613,301,247,433,283,829,689,991,249,406,472,926,1007,203,719,357,434,0,1023,212,424,519,337,442,526,614,163,81,373,223,540,457,716,428,326,192,162,312,235,444,446,291,568,594,403,164,920,642,345,432,141,12,384,794,324,826,113,954,470,822,377,987,381,535,71,266,624,936,676,786,295,616,328,547,817,778,772,329,179,634,353,773,282,74,24,650,257,934,565,673,137,804,629,103,226,574,885,299,429,937,621,330,243,89,951,323,251,313,558,881,142,1010,21,248,736,275,849,593,840,340,549,418,79,1008,720,86,145,167,582,501,611,941,533,45,521,751,147,566,358,348,756,899,195,88,523,258,53,366,148,242,48,111,788,800,3,846,845,901,618,538,834,798,274,873,585,783,746,150,206,892,452,85,636,425,747,955,87,22,347,411,851,774,730,293,149,1006,486,201,178,659,879,589,135,945,502,958,115,354,604,820,739,218,284,29,997,50,42,711,496,998,960,475,39,276,675,763,674,105,657,460,169,622,586,394,325,752,158,880,993,117,928,722,172,134,290,909,334,73,652,361,491,343,710,114,859,610,554,116,90,707,530,919,990,333,294,200,620,19,205,168,185,351,1000,797,775,304,390,451,176,189,534,923,5,805,106,771,221,563,296,367,484,214,96,156,222,401,553,82,577,216,6,653,669,733,667,1005,779,133,724,649,564,785,645,420,573,898,37,581,723,848,658,307,543,405,980,165,300,417,412,952,761,124,393,808,170,209,760,43,339,506,982,278,887,539,174,961,44,129,183,121,311,656,679,962,525,655,948,903,75,702,298,468,989,11,461,643,402,759,356,841,806,984,735,177,666,109,270,25,867,755,493,138,893,308,230,567,197,632,696,314,617,67,966,292,436,427,57,561,58,865,971,422,100,265,84,431,910,152,481,350,973,914,897,793,439,107,78,456,41,108,838,878,1014,322,836,904,210,705,802,680,409,458,338,464,732,62,660,360,277,253,139,781,992,320,316,584,737,839,963,607,234,261,833,891,932,748,344,310,268,889,69,154,795,572,157,289,146,469,792,944,211,388,471,76,175,713,908,309,228,54,695,161,708,864,596,229,232,31,180,382,902,160,548,931,815,386,957,630,155,956,77,542,400,490,728,194,38,868,410,27,336,688,370,271,191,80,977,193,571,734,527,245,97,690,269,600,391,40,352,623,378,601,556,20,823,812,10,662,587,512};

//S-box (S_5) (AES S-box)
int AES[]={99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22};

int nt, Nt;

int _tmain(int argc, _TCHAR* argv[])
{
	void fastwh(int *T, int *FW);

	int i,j,k,tmp,st;

	printf("\nSelect Sbox\n 1: S_1 (RSSB with nonlinearity 456)\n 2: S_2 (RSSB with nonlinearity 454)\n 3: S_3 (Concatenation with nonlinearity 456)\n 4: S_4 (Concatenation with nonlinearity 454)\n 5: S_5 (AES S-box)\n  ");
	scanf("%d",&st);

	if (st==1)
	{	for (i=0;i<N;i++)	S[i]=RSSB1[i];	Nt=N;	nt=n;	}
	else if (st==2)
	{	for (i=0;i<N;i++)	S[i]=RSSB2[i];	Nt=N;	nt=n;	}
	else if (st==3)
	{	for (i=0;i<N;i++)	S[i]=CONC1[i];	Nt=N;	nt=n;	}
	else if (st==4)
	{	for (i=0;i<N;i++)	S[i]=CONC2[i];	Nt=N;	nt=n;	}
	else if (st==5)
	{	for (i=0;i<256;i++)	S[i]=AES[i];	Nt=256;	nt=8;	}
	else
		return 0;


	//DDT: Difference distribution table.
	for (i=0;i<Nt;i++)
		for (j=0;j<Nt;j++)
			*(DDT+j*Nt+i)=0;

	for (i=0;i<Nt;i++)
		for (j=0;j<Nt;j++)
		{
			tmp=S[j^i]^S[j];
			*(DDT+tmp*Nt+i)=*(DDT+tmp*Nt+i)+1;
		}

	//DB: Binary values of all the inputs of an S-box.
	for (i=0;i<Nt;i++)
		for (j=0;j<nt;j++)
			*(DB+i*nt+j)=(i&(1<<(nt-1-j)))>>(nt-1-j);

	//Sb: Binary version of S.
	for (i=0;i<Nt;i++)
		for (j=0;j<nt;j++)
			*(Sb+i*nt+j)=*(DB+S[i]*nt+j);

	for (i=0;i<Nt;i++)
	{
		for (j=0;j<Nt;j++)
		{
			T[j]=0;
			for (k=0;k<nt;k++)
				T[j]=T[j]^((*(Sb+j*nt+k))*(*(DB+i*nt+k)));
		}
		//fastwh: Computes Walsh-Hadamard spectrum.
		fastwh(T,FW);
		for (j=0;j<Nt;j++)
			*(LAT+j*Nt+i)=FW[j]/2;
	}

	//Column and row numbers correspond to input and output sums, respectively.
	fprintf(outlat,"      ");
	for (i=0;i<Nt;i++)
		fprintf(outlat,"%4X ",i);
	fprintf(outlat,"\n");
	fprintf(outlat,"     ");
	for (i=0;i<Nt;i++)
		fprintf(outlat,"-----");
	for (i=0;i<Nt;i++)
	{
		fprintf(outlat,"\n%4X |",i);
		for (j=0;j<Nt;j++)
			fprintf(outlat,"%4d ",LAT[i*Nt+j]);
	}

	//Column and row numbers correspond to input and output differences, respectively.
	fprintf(outddt,"      ");
	for (i=0;i<Nt;i++)
		fprintf(outddt,"%4X ",i);
	fprintf(outddt,"\n");
	fprintf(outddt,"     ");
	for (i=0;i<Nt;i++)
		fprintf(outddt,"-----");
	for (i=0;i<Nt;i++)
	{
		fprintf(outddt,"\n%4X |",i);
		for (j=0;j<Nt;j++)
			fprintf(outddt,"%4d ",DDT[j*Nt+i]);
	}
	fclose(outlat);
	fclose(outddt);

	return 0;
}

void fastwh(int *T, int *FW)
{	
	int i,j,i1,i2,i3,k1=2,k2=Nt/2,k3=1,L1,temp1,temp2;
	for (i=0;i<Nt;i++)
		FW[i]=1-2*T[i];
	for (i1=0;i1<nt;i1++)  
	{
	   L1=1;
	   for (i2=0;i2<k2;i2++)
	   {
		  for (i3=0;i3<k3;i3++)
		  {
			 i=i3+L1-1; j=i+k3; 
		     temp1= FW[i]; temp2 = FW[j]; 
			 FW[i]=temp1+temp2;
		     FW[j]=temp1-temp2;
		  }
	      L1=L1+k1; 
	   }
	   k1=k1*2; k2=k2/2; k3=k3*2;
	}
}