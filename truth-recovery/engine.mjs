// engine.mjs -- pure frequentist NMA core EXTRACTED VERBATIM from index.html
// (LivingNMA). Stats helpers + matrix utils + buildNetwork + runNMACore copied
// unchanged; parseData (DOM) is NOT copied -- the harness feeds data directly.
// The SAME math the app ships is what the truth-recovery harness measures.

function normalCDF(z){const a1=0.254829592,a2=-0.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429,p=0.3275911;const s=z<0?-1:1;const x=Math.abs(z)/Math.sqrt(2);const t=1/(1+p*x);const y=1-(((((a5*t+a4)*t)+a3)*t+a2)*t+a1)*t*Math.exp(-x*x);return 0.5*(1+s*y)}
function normalQuantile(p){if(p<=0)return-Infinity;if(p>=1)return Infinity;if(p===0.5)return 0;const a=p<0.5?p:1-p;const t=Math.sqrt(-2*Math.log(a));const c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308;let z=t-(c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t);if(p<0.5)z=-z;return z}
function chiSqP(x,df){if(x<=0)return 1;return 1-regGammaP(df/2,x/2)}
function regGammaP(a,x){if(x<=0)return 0;const lga=lgamma(a);if(x<a+1){let s=1/a,t=1/a;for(let n=1;n<200;n++){t*=x/(a+n);s+=t;if(Math.abs(t)<1e-12*Math.abs(s))break}return s*Math.exp(-x+a*Math.log(x)-lga)}else{let f=1,b2=x+1-a,c2=1/1e-30,d2=1/b2;f=d2;for(let i=1;i<=200;i++){const an=-i*(i-a),bn=x+2*i+1-a;d2=bn+an*d2;if(Math.abs(d2)<1e-30)d2=1e-30;d2=1/d2;c2=bn+an/c2;if(Math.abs(c2)<1e-30)c2=1e-30;f*=d2*c2;if(Math.abs(d2*c2-1)<1e-10)break}return 1-f*Math.exp(-x+a*Math.log(x)-lga)}}
function lgamma(x){const g=7,c=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];if(x<0.5)return Math.log(Math.PI/Math.sin(Math.PI*x))-lgamma(1-x);x-=1;let a=c[0];const t=x+g+0.5;for(let i=1;i<g+2;i++)a+=c[i]/(x+i);return 0.5*Math.log(2*Math.PI)+(x+0.5)*Math.log(t)-t+Math.log(a)}
function tQuantile(p,df){
  if(df<=0)return normalQuantile(p);
  if(df>=300)return normalQuantile(p);
  const z=normalQuantile(p);
  const z2=z*z;
  // Hill's approximation for t quantile
  const g1=(z2+1)/4;
  const g2=(5*z2*z2+16*z2+3)/96;
  const g3=(3*z2*z2*z2+19*z2*z2+17*z2-15)/384;
  const g4=(79*z2*z2*z2*z2+776*z2*z2*z2+1482*z2*z2-1920*z2-945)/92160;
  return z+(g1/df)+(g2/(df*df))+(g3/(df*df*df))+(g4/(df*df*df*df));
}

// Matrix utilities
function matCreate(r,c,v=0){return Array.from({length:r},()=>Array(c).fill(v))}
function matMul(A,B){const r=A.length,c=B[0].length,k=B.length;const C=matCreate(r,c);for(let i=0;i<r;i++)for(let j=0;j<c;j++)for(let l=0;l<k;l++)C[i][j]+=A[i][l]*B[l][j];return C}
function matT(A){const r=A.length,c=A[0].length;const B=matCreate(c,r);for(let i=0;i<r;i++)for(let j=0;j<c;j++)B[j][i]=A[i][j];return B}
function matInv(M){
  const n=M.length;const A=M.map(r=>[...r]);const I=matCreate(n,n);
  for(let i=0;i<n;i++)I[i][i]=1;
  for(let i=0;i<n;i++){
    let maxR=i;for(let r=i+1;r<n;r++)if(Math.abs(A[r][i])>Math.abs(A[maxR][i]))maxR=r;
    [A[i],A[maxR]]=[A[maxR],A[i]];[I[i],I[maxR]]=[I[maxR],I[i]];
    const piv=A[i][i];if(Math.abs(piv)<1e-14)continue;
    for(let j=0;j<n;j++){A[i][j]/=piv;I[i][j]/=piv}
    for(let r=0;r<n;r++){if(r===i)continue;const f=A[r][i];for(let j=0;j<n;j++){A[r][j]-=f*A[i][j];I[r][j]-=f*I[i][j]}}
  }
  return I;
}
function buildNetwork(data){
  const trts=new Set();
  data.forEach(d=>{trts.add(d.t1);trts.add(d.t2)});
  const trtList=[...trts].sort();
  const trtIdx=Object.fromEntries(trtList.map((t,i)=>[t,i]));
  const T=trtList.length;

  // Adjacency
  const adj=matCreate(T,T,0);
  const edgeStudies={};
  for(const d of data){
    const i=trtIdx[d.t1],j=trtIdx[d.t2];
    adj[i][j]++;adj[j][i]++;
    const key=i<j?`${i}-${j}`:`${j}-${i}`;
    if(!edgeStudies[key])edgeStudies[key]=[];
    edgeStudies[key].push(d);
  }

  // Check connectivity (BFS)
  const visited=new Set([0]);
  const queue=[0];
  while(queue.length>0){
    const u=queue.shift();
    for(let v=0;v<T;v++){if(adj[u][v]>0&&!visited.has(v)){visited.add(v);queue.push(v)}}
  }
  const connected=visited.size===T;

  return{trtList,trtIdx,T,adj,edgeStudies,connected};
}

// Frequentist NMA using graph-theoretic approach (Rücker 2012)
function runNMACore(data,net,refIdx){
  const {trtList,trtIdx,T}=net;
  const K=data.length;// number of contrasts

  // Build design matrix X (K x T-1): each row encodes the contrast
  // We drop the reference treatment column
  const nonRefIdx=[];
  for(let i=0;i<T;i++)if(i!==refIdx)nonRefIdx.push(i);
  const X=matCreate(K,T-1);
  const Y=[];const W=[];

  for(let k=0;k<K;k++){
    const d=data[k];
    const i=trtIdx[d.t1],j=trtIdx[d.t2];
    // Convention: y = effect of t1 vs t2
    Y.push(d.y);
    W.push(1/d.v);
    // X row: +1 for t1's column (if not ref), -1 for t2's column (if not ref)
    for(let c=0;c<T-1;c++){
      const tIdx=nonRefIdx[c];
      if(tIdx===i)X[k][c]=1;
      else if(tIdx===j)X[k][c]=-1;
      else X[k][c]=0;
    }
  }

  // Weighted least squares: beta = (X'WX)^-1 X'WY
  const Wmat=matCreate(K,K,0);
  for(let k=0;k<K;k++)Wmat[k][k]=W[k];
  const Xt=matT(X);
  const XtW=matMul(Xt,Wmat);
  const XtWX=matMul(XtW,X);
  const XtWXinv=matInv(XtWX);
  const Yvec=Y.map(y=>[y]);
  const XtWY=matMul(XtW,Yvec);
  const beta=matMul(XtWXinv,XtWY).map(r=>r[0]);

  // Variance-covariance of beta = (X'WX)^-1
  const vcov=XtWXinv;

  // Q statistic
  const fitted=matMul(X,beta.map(b=>[b])).map(r=>r[0]);
  let Q=0;
  for(let k=0;k<K;k++){Q+=W[k]*(Y[k]-fitted[k])*(Y[k]-fitted[k])}
  const df=K-(T-1);
  const Qp=df>0?chiSqP(Q,df):1;

  // DL heterogeneity
  let tau2=0;
  if(df>0&&Q>df){
    let sumW=0,trWXtWXinvW=0;
    for(let k=0;k<K;k++)sumW+=W[k];
    // trace(W - W*X*(X'WX)^-1*X'W) = sum(W) - trace(XtWX_inv * XtW * W * X)
    // Simplified: C = sumW - trace(...)
    const WX=matMul(Wmat,X);
    const H=matMul(matMul(X,XtWXinv),Xt);// Hat matrix without W
    let trHW=0;
    for(let k=0;k<K;k++)trHW+=H[k][k]*W[k]*W[k]/W[k];// simplified
    const C=sumW-trHW;
    tau2=C>0?Math.max(0,(Q-df)/C):0;
  }

  // Build all pairwise results
  const comparisons=[];
  const alpha=0.05;
  const zc=normalQuantile(1-alpha/2);

  // vs reference comparisons
  for(let c=0;c<T-1;c++){
    const trt=trtList[nonRefIdx[c]];
    const est=beta[c];
    const se=Math.sqrt(vcov[c][c]);
    const z=est/se;
    const p=2*(1-normalCDF(Math.abs(z)));
    comparisons.push({t1:trt,t2:trtList[refIdx],est,se,ci:[est-zc*se,est+zc*se],p,z});
  }

  // All pairwise (for league table)
  const allPairs=matCreate(T,T,null);
  for(let i=0;i<T;i++){
    for(let j=0;j<T;j++){
      if(i===j)continue;
      // Effect of trtList[i] vs trtList[j]
      if(i===refIdx){
        // -beta[j's index in nonRef]
        const jc=nonRefIdx.indexOf(j);
        if(jc>=0){
          const est=-beta[jc];const se=Math.sqrt(vcov[jc][jc]);
          allPairs[i][j]={est,se,ci:[est-zc*se,est+zc*se],p:2*(1-normalCDF(Math.abs(est/se)))};
        }
      }else if(j===refIdx){
        const ic=nonRefIdx.indexOf(i);
        if(ic>=0){
          const est=beta[ic];const se=Math.sqrt(vcov[ic][ic]);
          allPairs[i][j]={est,se,ci:[est-zc*se,est+zc*se],p:2*(1-normalCDF(Math.abs(est/se)))};
        }
      }else{
        const ic=nonRefIdx.indexOf(i),jc=nonRefIdx.indexOf(j);
        if(ic>=0&&jc>=0){
          const est=beta[ic]-beta[jc];
          const se=Math.sqrt(vcov[ic][ic]+vcov[jc][jc]-2*vcov[ic][jc]);
          allPairs[i][j]={est,se,ci:[est-zc*se,est+zc*se],p:2*(1-normalCDF(Math.abs(est/se)))};
        }
      }
    }
  }

  // P-scores
  const pscores=[];
  for(let i=0;i<T;i++){
    let score=0,count=0;
    for(let j=0;j<T;j++){
      if(i===j||!allPairs[i][j])continue;
      const z2=allPairs[i][j].est/allPairs[i][j].se;
      // P(trt i better than trt j) — assumes lower is better for log scale
      score+=normalCDF(-z2);// negative z means i is better (lower estimate)
      count++;
    }
    pscores.push({trt:trtList[i],idx:i,pscore:count>0?score/count:0.5});
  }
  pscores.sort((a,b)=>b.pscore-a.pscore);

  // Node-splitting (local inconsistency)
  const nodesplit=[];
  const directPairs=new Set();
  data.forEach(d=>{
    const i=trtIdx[d.t1],j=trtIdx[d.t2];
    const key=Math.min(i,j)+'-'+Math.max(i,j);
    directPairs.add(key);
  });

  for(const key of directPairs){
    const [i,j]=key.split('-').map(Number);
    if(!allPairs[i][j])continue;
    // Direct evidence
    const directStudies=data.filter(d=>(trtIdx[d.t1]===i&&trtIdx[d.t2]===j)||(trtIdx[d.t1]===j&&trtIdx[d.t2]===i));
    if(directStudies.length===0)continue;
    let dSumW=0,dSumWY=0;
    for(const ds of directStudies){
      const w=1/ds.v;
      const sign=(trtIdx[ds.t1]===i)?1:-1;
      dSumW+=w;dSumWY+=w*ds.y*sign;
    }
    const directEst=dSumWY/dSumW;
    const directSE=1/Math.sqrt(dSumW);
    // Indirect = NMA estimate (consistency model)
    const nmaEst=allPairs[i][j].est;
    const nmaSE=allPairs[i][j].se;
    // Indirect estimate from NMA minus direct
    // Approximate: indirect = (nmaEst/nmaSE^2 - directEst/directSE^2) / (1/nmaSE^2 - 1/directSE^2)
    const wNMA=1/(nmaSE*nmaSE),wDir=1/(directSE*directSE);
    if(wNMA<=wDir)continue;// can't separate
    const indirectEst=(nmaEst*wNMA-directEst*wDir)/(wNMA-wDir);
    const indirectSE=1/Math.sqrt(Math.abs(wNMA-wDir));

    const diff=directEst-indirectEst;
    const diffSE=Math.sqrt(directSE*directSE+indirectSE*indirectSE);
    const pVal=2*(1-normalCDF(Math.abs(diff/diffSE)));

    nodesplit.push({
      t1:trtList[i],t2:trtList[j],
      direct:{est:directEst,se:directSE},
      indirect:{est:indirectEst,se:indirectSE},
      diff,diffSE,p:pVal
    });
  }

  // === Prediction intervals for all pairs ===
  // PI = est +/- t_{k_designs - T + 1, alpha/2} * sqrt(se^2 + tau2)
  // k_designs = number of unique direct comparisons (designs)
  const kDesigns=Object.keys(net.edgeStudies).length;
  const piDf=Math.max(1, kDesigns - T + 1);
  const tCrit=tQuantile(1-alpha/2, piDf);
  const allPI=matCreate(T,T,null);
  for(let i=0;i<T;i++){
    for(let j=0;j<T;j++){
      if(i===j||!allPairs[i][j])continue;
      const p=allPairs[i][j];
      const piSE=Math.sqrt(p.se*p.se+tau2);
      allPI[i][j]={lo:p.est-tCrit*piSE, hi:p.est+tCrit*piSE};
    }
  }

  // === Hat matrix H = X (X'WX)^-1 X'W for contribution matrix ===
  const Hmat=matMul(matMul(X,XtWXinv),XtW);

  // === Store design matrix and weights for advanced analyses ===
  const Xmat=X, Yvec2=Y, Wvec=W, fittedVals=fitted;

  return{beta,vcov,Q,Qp,df,tau2,comparisons,allPairs,pscores,nodesplit,nonRefIdx,
    allPI,piDf,tCrit,Hmat,Xmat,Yvec:Yvec2,Wvec,fittedVals,kDesigns};
}


export { buildNetwork, runNMACore, matCreate, matMul, matT, matInv, normalCDF, normalQuantile, chiSqP };
