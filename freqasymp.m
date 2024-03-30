function [A, F, lgw1, lgws]=freqasymp(sys, lgw, flg)
%freqasymp Вычисление асимптотических логарифмических АФЧХ
%sys-lti модель SISO типа zpk или tf; после каждого комплексного
%нуля или полюса должно сразу следовать сопряженное ему значение
%lgw-вектор десятичных логарифмов круговых частот или
%границы диапазона: [log10(wmin) log10(wmax)]
%flg-если 1, то вычисляется асимптотическая ФЧХ
%если 0 или отсутствует, то просто ФЧХ
%A - АЧХ (дБ)
%F - ФЧХ (градусы)
%lgwl - логарифмы частот, на которых выч. A и F
%lgws - логарифмы сопрягающих частот
if ~((isa(sys,'zpk')||isa(sys,'tf'))&&issiso(sys))
    error('Требуется модель SISO типа "zpk" или "tf".'),end
if isa(sys,'tf'), sys = zpk(sys); end
lgw = lgw(:); ws=[]; lgws=[]; kn=2*log(10); kw=pi/kn;
if (nargin == 3 && ~flg)||nargin==2, flag=1; else flag=0; end
if flag, kdf = 5; df = (lgw(end) - lgw(1))/kdf; end
zp = [sys.z{:}; sys.p{:}];
if ~isempty(zp)
    jr = find(~imag(zp));
    if ~isempty(jr)
        zpr = zp(jr);
        jr1 = find(zpr);
        if ~isempty(jr1)
            zpr1 = zpr(jr1);
            w = log10(abs(zpr1));
            if flag, df = min(df, kw/kdf);
        end
            lgws = [lgws; w];
            ws = [ws; w; w-kw; w+kw];
    end
end
ji = find(imag(zp));
if ~isempty(ji)
    zpi = zp(ji); zpi = zpi(1:2:end);
    w = log10(abs(zpi));
    dw = kw*abs(real(zpi))./abs(zpi);
        if flag, jj = find(dw);
            if ~isempty(jj), df = min(df,min(dw(jj))/kdf); end
end
lgws = [lgws; w];
ws = [ws; w; w-dw; w+dw];
end
end
if flag
lgw1 = sort([lgw; ws; (lgw(1):df:lgw(end))']);
ww = 10.^lgw1;
else
lgw1 = sort([lgw; ws]);
end
k = sys.k;
A = log10(abs(k)) + zeros(length(lgw1),1);
F = atan2(0, k) + zeros(length(lgw1),1);
p = sys.p{:}; % плюсы
if ~isempty(p)
jr = find(~imag(p));
if ~isempty(jr)
pr = p(jr);
j0 = find(~pr);
if ~isempty(j0)
A = A - lgw1*length(j0);
F = F - pi/2*length(j0);
end
jr1 = find(pr);
if ~isempty(jr1)
pr1 = pr(jr1);
w = log10(abs(pr1));
for j = 1:length(pr1)
jw = find(lgw1 <= w(j));
A(jw) = A(jw)-w(j); jn = jw(end)+1;
A(jn:end) = A(jn:end)- lgw1(jn:end);
if flag, F = F - atan2(ww, -pr1(j));
else
jw = find(lgw1 > w(j)-kw&lgw1<= w(j)+kw);
if pr1(j)<0
F(jw)=F(jw)-(kn*(lgw1(jw)-w(j))+pi)*0.25;
else
jn = jw(1)-1;
F(1:jn) = F(1:jn) - pi;
F(jw)=F(jw)+(kn*(lgw1(jw)-w(j))-3*pi)*0.25;
end
jn = jw(end)+1;
F(jn:end) = F(jn:end) - pi*0.5;
end
end
end
end
ji = find(imag(p));
if ~isempty(ji)
pm = p(ji); pm = pm(1:2:end);
w = log10(abs(pm));
dw = kw*abs(real(pm))./abs(pm);
for j = 1:length(pm)
jw = find(lgw1 <= w(j));
A(jw) = A(jw) - 2*w(j); jn =jw(end)+1;
A(jn:end) = A(jn:end) - 2*lgw1(jn:end);
if flag
if real(pm(j))
F= F-atan2(ww-imag(pm(j)),-real(pm(j)));
F= F-atan2(ww+imag(pm(j)),-real(pm(j)));
end
else
jw=find(lgw1 > w(j)-dw(j)& lgw1<= w(j)+dw(j));
jn = jw(end) + 1;
if real(pm(j))<0
F(jw)=F(jw)-((lgw1(jw)-w(j))/dw(j)+1)*pi*0.5;
F(jn:end) = F(jn:end) - pi;
elseif real(pm(j))>0
F(jw)=F(jw)+((lgw1(jw)- w(j))/dw(j)+1)*pi*0.5;
F(jn:end) = F(jn:end) + pi;
end
end
if real(pm(j))== 0
jw = find(lgw1 > w(j)); F(jw) = F(jw) - pi;
end
end
end
end
z = sys.z{:}; % нули
if ~isempty(z)
jr = find(~imag(z));
if ~isempty(jr)
zr = z(jr);
j0 = find(~zr);
if ~isempty(j0)
A = A + lgw1*length(j0);
F = F + pi/2*length(j0);
end
jr1 = find(zr);
if ~isempty(jr1)
    zr1 = zr(jr1);
w = log10(abs(zr1));
for j = 1:length(zr1)
jw = find(lgw1 <= w(j));
A(jw) = A(jw) + w(j); jn = jw(end) + 1;
A(jn:end) = A(jn:end) + lgw1(jn:end);
if flag, F = F + atan2(ww, -zr1(j));
else
jw = find(lgw1 > w(j)-kw & lgw1 <= w(j)+kw);
if zr1(j) < 0
F(jw)=F(jw)+(kn*(lgw1(jw)-w(j)) + pi)*0.25;
else
jn = jw(1) - 1;
F(1:jn) = F(1:jn) + pi;
F(jw)=F(jw)-(kn*(lgw1(jw)-w(j))-3*pi)*0.25;
end
jn = jw(end) + 1;
F(jn:end) = F(jn:end) + pi*0.5;
end
end
end
end
ji = find(imag(z));
if ~isempty(ji)
zm = z(ji); zm = zm(1:2:end);
w = log10(abs(zm));
dw = kw*abs(real(zm))./abs(zm);
for j = 1:length(zm)
jw = find(lgw1 <= w(j));
A(jw) = A(jw) + 2*w(j); jn = jw(end) + 1;
A(jn:end) = A(jn:end) + 2*lgw1(jn:end);
if flag
if real(zm(j))
F = F + atan2(ww-imag(zm(j)),-real(zm(j)));
F = F + atan2(ww+imag(zm(j)),-real(zm(j)));
end
else
jw=find(lgw1 >w(j)-dw(j) & lgw1<= w(j)+dw(j));
jn = jw(end) + 1;
if real(zm(j))<0
F(jw)=F(jw)+((lgw1(jw)-w(j))/dw(j)+1)*pi*0.5;
F(jn:end) = F(jn:end) + pi;
elseif real(zm(j))>0
F(jw)=F(jw)-((lgw1(jw)-w(j))/dw(j)+1)*pi*0.5;
F(jn:end) = F(jn:end)-pi;
end
end
if real(zm(j))==0
jw = find(lgw1 > w(j)); F(jw) = F(jw) + pi;
end
end
end
end
j = find(lgw1>=lgw(1) & lgw1<=lgw(end));
lgw1 = lgw1(j); A = A(j); F = F(j);
if flag, F = unwrap(F); end
A = 20*A; F = F*180/pi;