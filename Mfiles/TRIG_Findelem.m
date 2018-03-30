function TheReply=TRIG_Findelem(ListTrig, ListNode, PT)

nbTrig=size(ListTrig, 1);
LSi=[1 2 3];
LSj=[2 3 1];
LSk=[3 1 2];

PTmat=zeros(nbTrig,2);
PTmat(:,1)=PT(1,1);
PTmat(:,2)=PT(1,2);



CResult=zeros(nbTrig, 3);
for iPos=1:3
  iNode=ListTrig(:,LSi(1,iPos));
  jNode=ListTrig(:,LSj(1,iPos));
  kNode=ListTrig(:,LSk(1,iPos));
  V1=ListNode(jNode,:)-ListNode(iNode,:);
  V2=ListNode(kNode,:)-ListNode(iNode,:);
  Vdiff=PTmat-ListNode(iNode,:);
  det1=V1(:,1).*V2(:,2)-V1(:,2).*V2(:,1);
  det2=V1(:,1).*Vdiff(:,2)-V1(:,2).*Vdiff(:,1);
  ScalProd=det1.*det2 > 0;
  CResult(:,iPos)=ScalProd;
end;
CBelong=CResult(:,1).*CResult(:,2).*CResult(:,3);
ListReply=find(CBelong == 1);
nbBelong=size(ListReply,1);
if (nbBelong > 1)
  disp(['The vertex belongs to ' num2str(nbBelong) ' triangles']);
end;
if (nbBelong == 0)
  TheReply=-1;
else
  TheReply=ListReply(1,1);
end;

