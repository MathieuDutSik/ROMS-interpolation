function ListTrigNr=TRIG_Findelems(ListTrig, ListNode, ListPT)

nbNode=size(ListNode,1);
nbPoint=size(ListPT,1);
nbTrig=size(ListTrig,1);
ListNodeLON=ListNode(:,1);
ListNodeLAT=ListNode(:,2);


%method=1; % griddata
%method=2; % Rich signell nearest
method=3; % Rich signell nearest + preselection of good nodes
if (method == 1)
  disp('Using method 1');
  iListNode=1:nbNode;
  iListPT=griddata(ListNodeLON, ListNodeLAT, iListNode, ...
		   ListPT(:,1), ListPT(:,2), 'nearest');
elseif (method == 2)
  disp('Using method 2');
  iListPT=zeros(nbPoint,1);
  for ipt=1:nbPoint
    idx=nearxy(ListNodeLON, ListNodeLAT, ...
	       ListPT(ipt, 1), ListPT(ipt,2));
    iListPT(ipt,1)=idx;
  end;
elseif (method == 3)
  disp('Using method 3');
  deltaLonLat=0.5;
  MinLon=min(ListPT(:, 1))-deltaLonLat;
  MaxLon=max(ListPT(:, 1))+deltaLonLat;
  MinLat=min(ListPT(:, 2))-deltaLonLat;
  MaxLat=max(ListPT(:, 2))+deltaLonLat;
  K1=find(ListNodeLON > MinLon);
  K2=find(ListNodeLON < MaxLon);
  K3=find(ListNodeLAT > MinLat);
  K4=find(ListNodeLAT < MaxLat);
  K12=intersect(K1, K2);
  K34=intersect(K3, K4);
  K1234=intersect(K12, K34);
  Lpos=zeros(nbNode, 1);
  for iNode=1:nbNode
    Lpos(iNode, 1)=iNode;
  end;
  nbNodeRed=size(K1234, 1);
  disp(['Reduce from nbNode=' num2str(nbNode) ...
	' to nbNodeRed=' num2str(nbNodeRed)]);
  LposRed=Lpos(K1234);
  iListPT=zeros(nbPoint,1);
  ListNodeLONred=ListNodeLON(K1234);
  ListNodeLATred=ListNodeLAT(K1234);
  for ipt=1:nbPoint
    idxred=nearxy(ListNodeLONred, ListNodeLATred, ...
	       ListPT(ipt, 1), ListPT(ipt,2));
    idx=LposRed(idxred);
    iListPT(ipt,1)=idx;
  end;
else
  disp('You did not choosed method correctly, please correct');
  keyboard;
end;
disp('nearest interpolation done');
disp(['nbTrig=' num2str(nbTrig)]);
[A,B]=meshgrid(1:nbTrig,1:3);
Ap=A';
Bp=B';

ListTrigNr=zeros(nbPoint,1);
nbCase1=0;
nbCase2=0;
nbCase3=0;
%Adebug=zeros(3,3);
for iPoint=1:nbPoint
  PT=ListPT(iPoint,:);
  iNodeRel=iListPT(iPoint,1);
  K=find(ListTrig == iNodeRel);
  ListSel=Ap(K);
  ListTrigSel=ListTrig(ListSel, :);
  TheReply=TRIG_Findelem(ListTrigSel, ListNode, PT);
  if (TheReply == -1)
    TheReply=TRIG_Findelem(ListTrig, ListNode, PT);
    if (TheReply == -1)
      nbCase3=nbCase3+1;
    else
      nbCase2=nbCase2+1;
    end;
    RightTrig=TheReply;
  else
    nbCase1=nbCase1+1;
    RightTrig=ListSel(TheReply);
  end;
%  if (RightTrig > 0)
%    iNode1=ListTrig(RightTrig,1);
%    iNode2=ListTrig(RightTrig,2);
%    iNode3=ListTrig(RightTrig,3);
%    Adebug(1,:)=[1 ListNode(iNode1,1) ListNode(iNode1,2)];
%    Adebug(2,:)=[1 ListNode(iNode2,1) ListNode(iNode2,2)];
%    Adebug(3,:)=[1 ListNode(iNode3,1) ListNode(iNode3,2)];
%    Bdebug=[1 PT(1,1) PT(1,2)];
%    eSol=Bdebug*inv(Adebug);
%    if (eSol(1,1) < -0.05 || eSol(1,1) > 1.05 || ...
%	eSol(1,2) < -0.05 || eSol(1,2) > 1.05 || ...
%	eSol(1,3) < -0.05 || eSol(1,3) > 1.05)
%      disp(['Please debug here in TRIG_Findelems.m']);
%      keyboard;
%    end;
%  end;
  ListTrigNr(iPoint,1)=RightTrig;
end;
disp(['ansatz success : ' num2str(nbCase1) '  cases']);
disp(['direct success : ' num2str(nbCase2) '  cases']);
disp(['      failures : ' num2str(nbCase3) '  cases']);
