
#ifndef SYSTEM_H_
#define SYSTEM_H_

#ifdef _WIN32
#define cpuid(info,x)    __cpuidex(info,x,0)
#else
void cpuid(int CPUInfo[4],int InfoType) {
    __asm__ __volatile__ (
        "cpuid":
        "=a" (CPUInfo[0]),
        "=b" (CPUInfo[1]),
        "=c" (CPUInfo[2]),
        "=d" (CPUInfo[3]) :
        "a" (InfoType), "c" (0)
    );
}
#endif

bool check_SSSE3()
{
    int info[4];
    cpuid(info, 0);
    int nids = info[0];
    if (nids >= 1) {
    	cpuid(info, 1);
    	return (info[2] & (1<<9)) != 0;
    }
    return false;
}

#endif /* SYSTEM_H_ */
