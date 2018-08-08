# @Date: 11 January 2018
# @Modified: -
# @Author: Myrthe Jager
# @Description: Determine whether SNV data is normally distributed
# Abbreviations: -




#1 Manually enter the mutation numbers (from Blokzijl, et al. Nature 2016)
mutations_si <- c(246, 245, 304,591, 506, 614, 1461, 1838, 1571, 2497, 2613, 3516, 3512, 1957)
mutations_liver <- c(771,888, 1292, 1351, 1495, 1273, 1845, 1504, 1577, 1919)

#2 Manually enter the age of the donors (from Blokzijl, et al. Nature 2016)
age_si <- c(3,3,3,8,8,8,44,45,45,70,74,78,87,87)
age_liver <- c(30,30,41,41,46,46,53,53,53,55)

#3 Calculate mutation rate
rate.liver <- mutations_liver/age_liver
rate.si <- mutations_si/age_si

#4 Test for normality
shapiro.test(rate.liver)
shapiro.test(rate.si)