clear; close all; clc

load(".\TimoShear\strainSection1ele1.out")
load(".\TimoShear\stressSection1ele1.out")
load(".\TimoShear\strainSection2ele1.out")
load(".\TimoShear\stressSection2ele1.out")
load(".\TimoShear\strainSection3ele1.out")
load(".\TimoShear\stressSection3ele1.out")
load(".\TimoShear\strainSection4ele1.out")
load(".\TimoShear\stressSection4ele1.out")
load(".\TimoShear\strainSection5ele1.out")
load(".\TimoShear\stressSection5ele1.out")

figure
for i = 1:5
        subplot(2,5,i)
        plot(strainSection1ele1(:,2*i-1), stressSection1ele1(:,2*i-1),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," normal stress");
        title(str)

        subplot(2,5,5+i)
        plot(strainSection1ele1(:,2*i), stressSection1ele1(:,2*i),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," shear stress");
        title(str)
end

figure
for i = 1:5
        subplot(2,5,i)
        plot(strainSection2ele1(:,2*i-1), stressSection2ele1(:,2*i-1),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," normal stress");
        title(str)

        subplot(2,5,5+i)
        plot(strainSection2ele1(:,2*i), stressSection2ele1(:,2*i),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," shear stress");
        title(str)
end

figure
for i = 1:5
        subplot(2,5,i)
        plot(strainSection3ele1(:,2*i-1), stressSection3ele1(:,2*i-1),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," normal stress");
        title(str)

        subplot(2,5,5+i)
        plot(strainSection3ele1(:,2*i), stressSection3ele1(:,2*i),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," shear stress");
        title(str)
end

figure
for i = 1:5
        subplot(2,5,i)
        plot(strainSection4ele1(:,2*i-1), stressSection4ele1(:,2*i-1),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," normal stress");
        title(str)

        subplot(2,5,5+i)
        plot(strainSection4ele1(:,2*i), stressSection4ele1(:,2*i),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," shear stress");
        title(str)
end

figure
for i = 1:5
        subplot(2,5,i)
        plot(strainSection5ele1(:,2*i-1), stressSection5ele1(:,2*i-1),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," normal stress");
        title(str)

        subplot(2,5,5+i)
        plot(strainSection5ele1(:,2*i), stressSection5ele1(:,2*i),'r');
        grid on
        xlabel('strain')
        ylabel('stress')
        str = strcat("fiber ",num2str(i)," shear stress");
        title(str)
end
