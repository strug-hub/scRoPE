// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <Rmath.h>
// #include <unsupported/Eigen/SpecialFunctions>


// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::COLAMDOrdering<int> CO;

//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double pmg_ll_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                    const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                    const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                    const Eigen::VectorXi & posindy_c, const int nind_c, const int k_c,
                    const Eigen::VectorXd & beta_c,const double sigma_c)
{

    double term1 = 0;

    Eigen::VectorXd xtb = offset_c + X_c*beta_c;

    Eigen::VectorXd extb = xtb.array().exp();

    Eigen::VectorXd cumsumxtb(k_c);
    for(int i=0;i<k_c;i++)
    {cumsumxtb(i)= extb.segment(fid_c(i),fid_c(i+1)-fid_c(i)).sum();}

    double exps = exp(sigma_c);
    double alpha = 1/(exps-1);
    double lambda = 1/(sqrt(exps)*(exps-1));

    unsigned int nelem = posind_c.size();
    Eigen::VectorXd cumsumy_ca = cumsumy_c.array() + alpha;
    for(unsigned int i=0;i<nelem;i++)
    {term1 += lgamma(cumsumy_ca(posind_c(i)));}
    term1 -= nelem*lgamma(alpha);

    nelem = posindy_c.size();
    for(unsigned int i=0;i<nelem;i++)
    {term1 += xtb(posindy_c(i))*Y_c(i);}

    term1 += k_c*alpha*log(lambda);

    Eigen::ArrayXd cumsumxtb_t = cumsumxtb.array()+lambda;
    term1 -= cumsumy_ca.dot(cumsumxtb_t.log().matrix());

    term1 = term1*(-1);
    return term1;

}


//
// [[Rcpp::export]]
Eigen::VectorXd pmg_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                              const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                              const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                              const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                              const Eigen::VectorXd & beta_c,const double sigma_c)
{
    Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+1);

    Eigen::VectorXd xtb = offset_c + X_c*beta_c;
    Eigen::VectorXd extb = xtb.array().exp();
    Eigen::ArrayXd cumsumxtb(k_c);
    Eigen::ArrayXd len(k_c);
    for(int i=0;i<k_c;i++)
    {
        len(i) = fid_c[i+1]-fid_c[i];
        cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
    }

    double exps = exp(sigma_c);
    double alpha = 1/(exps-1);
    double lambda = 1/(sqrt(exps)*(exps-1));

    double alpha_pr = -exps/pow(exps-1,2);
    double lambda_pr = (1-3*exps)/(2*sqrt(exps)*pow(exps-1,2));

    Eigen::ArrayXd ystar = cumsumy_c.array() + alpha;
    Eigen::ArrayXd mustar = cumsumxtb + lambda;
    Eigen::ArrayXd ymustar = ystar/mustar;

    Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
    Eigen::MatrixXd xexb_f(k_c,nb_c);
    for(int i=0;i<k_c;i++)
    {
        xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }

    Eigen::VectorXd db = Eigen::VectorXd::Zero(nb_c);
    db -= xexb_f.transpose()*ymustar.matrix();
    unsigned int nelem = posindy_c.size();
    for(unsigned int i=0;i<nelem;i++)
    {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

    double dtau = 0;
    double ldm = log(lambda)*k_c - mustar.log().sum();
    double adlmy = alpha/lambda*k_c - ymustar.sum();
    dtau += alpha_pr*ldm + lambda_pr*adlmy;

    der.segment(0,nb_c) = db;
    der[nb_c] = dtau;

    return (-1)*der;
}

// [[Rcpp::export]]
Eigen::MatrixXd pmg_hes_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::VectorXd & offset_c,
                              const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                              const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                              const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                              const Eigen::VectorXd & beta_c,const double sigma_c)
{
  
  Eigen::VectorXd xtb = offset_c + X_c*beta_c;
  Eigen::VectorXd extb = xtb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }
  
  double exps = exp(sigma_c);
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double log_lambda = log(lambda);
  
  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);
  double alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;
  double lambda_dpr = -3*exps/(2*exps_s*exps_m) + (3*exps-1)/(4*exps_s*exps_m) + (3*exps-1)*exps_s/(exps_m*(exps-1));
  
  Eigen::ArrayXd ystar = cumsumy_c.array() + alpha;
  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::ArrayXd ymustar = ystar/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::ArrayXd imustar = 1/mustar;
  
  Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
  }
  
  Eigen::MatrixXd hes(nb_c+1,nb_c+1); 
  Eigen::VectorXd tempb;
  for(int i=0;i<nb_c;i++)
  {
    for(int j=i;j<nb_c;j++)
    {
        tempb = xexb.col(i).cwiseProduct(X_c.col(j));
        double tempc = 0;
        for(int k=0;k<k_c;k++)
        {
          tempc += ymumustar(k)*xexb_f(k,i)*xexb_f(k,j) - ymustar(k)*tempb.segment(fid_c(k),len(k)).sum();
        }
        hes(i,j) = tempc;
        if(i!=j)
        {hes(j,i) = tempc;}
    }
  }
  tempb = lambda_pr*ymumustar-alpha_pr*imustar;
  hes.block(0,nb_c,nb_c,1) = xexb_f.transpose()*tempb;
  hes.block(nb_c,0,1,nb_c) = hes.block(0,nb_c,nb_c,1).transpose();
   
  double hes_sigma = alpha_dpr*log_lambda + 2*alpha_pr*lambda_pr/lambda - alpha*pow(lambda_pr/lambda,2) + alpha/lambda*lambda_dpr;
  hes_sigma = hes_sigma*k_c;
  hes_sigma = hes_sigma - (2*alpha_pr*lambda_pr*imustar.sum() + alpha_dpr*mustar.log().sum() + lambda_dpr*ymustar.sum() - lambda_pr*lambda_pr*ymumustar.sum());
  hes(nb_c,nb_c) = hes_sigma;
  
  return hes;
}


//
// [[Rcpp::export]]
double ptmg_ll_eigen(const Eigen::MatrixXd & X_c, const Eigen::VectorXd & offset_c,
                     const Eigen::VectorXd & Y_c,
                     const Eigen::VectorXd & fam_c, const Eigen::VectorXi & fid_c,
                     const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                     const Eigen::VectorXi & posindy_c, const int nind_c, const int k_c,
                     const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double term1 = 0;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();

  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}

  extb = extb.array().exp();

  Eigen::VectorXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  // unsigned int nelem_1 = posind_c.size();
  Eigen::VectorXd cumsumy_ca = cumsumy_c.array() + alpha;

  Eigen::ArrayXd cumsumxtb_t = cumsumxtb.array() + lambda;
  term1 -= cumsumy_ca.dot(cumsumxtb_t.log().matrix());

  term1 += k_c*alpha*log(lambda);
  term1 += nind_c*gamma*log(gamma);

  Eigen::VectorXd sum_elgcp(nind_c);
  cumsumy_ca = (cumsumy_ca.array()/cumsumxtb_t).matrix();


  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = cumsumy_ca(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += cumsumy_ca.dot(cumsumxtb.matrix());
  // term1 += gamma/phi*sum_elgcp.sum();

  // sum_elgcp = (phi+sum_elgcp.array()).log();
  sum_elgcp = (sum_elgcp.array()+gamma).log();
  term1 -= gamma*sum_elgcp.sum();

  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}

  // term1 = term1*(-1);

  return term1;

}


//
// [[Rcpp::export]]
Eigen::VectorXd ptmg_der_eigen(const Eigen::MatrixXd & X_c, const Eigen::VectorXd & offset_c,
                               const Eigen::VectorXd & Y_c, const Eigen::VectorXd & fam_c, const Eigen::VectorXi & fid_c,
                               const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                               const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                               const Eigen::VectorXd & beta_c,const Eigen::VectorXd & sigma_c) {

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::VectorXd xtb = offset_c + X_c*beta_c;
  Eigen::VectorXd extb = xtb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  for(int i=0;i<k_c;i++)
  {cumsumxtb(i)= extb.segment(fid_c(i),fid_c(i+1)-fid_c(i)).sum();}

  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double alpha_pr = -exps/pow(exps-1,2);
  double lambda_pr = (1-3*exps)/(2*sqrt(exps)*pow(exps-1,2));

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;
  Eigen::ArrayXd mustar = lambda + cumsumxtb;
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;

  // Eigen::VectorXd gstar = Y_c.array() + gamma;
  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}
  Eigen::VectorXd gstar_phiymustar(nind_c);
  Eigen::ArrayXd len(k_c);
  double slpey = 0;

  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    gstar_phiymustar.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }
  gstar_phiymustar = gamma + gstar_phiymustar.array();
  slpey = gstar_phiymustar.array().log().sum();
  gstar_phiymustar = gstar.array()/gstar_phiymustar.array();

  Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    // len(i) = fid_c[i+1]-fid_c[i];
    xexb_f.row(i) = xexb.block(fid_c[i],0,len(i),nb_c).colwise().sum();
  }
  xexb_f.transposeInPlace();

  Eigen::MatrixXd dbeta_41(k_c,nb_c);
  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
  }
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd ymumustar_csxtb = ymumustar*cumsumxtb;
  Eigen::ArrayXd ymumustar_dbeta = ymumustar*dbeta_42;
  // Eigen::VectorXd db = - pow(sexp2,2)*(xexb_f*ymumustar_csxtb.matrix()) - dbeta_41*ymustar.matrix() + sexp2*(xexb_f*(ymumustar*dbeta_42).matrix());
  Eigen::VectorXd db = xexb_f*(ymumustar_dbeta - ymumustar_csxtb).matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  double dtau = 0;
  double ldm = log(lambda)*k_c - mustar.log().sum();
  double adlmy = alpha/lambda*k_c - ymustar.sum();
  dtau += alpha_pr*(cumsumxtb/mustar).sum() - lambda_pr*ymumustar_csxtb.sum();
  dtau -= alpha_pr*(dbeta_42/mustar).sum() - lambda_pr*ymumustar_dbeta.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  double dtau2 = log(gamma)*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return (-1)*der;
}


// [[Rcpp::export]]
Eigen::MatrixXd call_cumsumy(const Eigen::MappedSparseMatrix<double> count, const Eigen::VectorXi & fid, const int k, const int ng)
{

  Eigen::MatrixXd cumsumy(ng,k);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(ng);
  int temp_k = 0;

  for (int i=0; i<count.outerSize(); ++i)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(count,i); it; ++it)
    {
      if(it.value()>0)
      {
        temp[it.row()] += it.value();
      }
    }
    if(i==(fid[temp_k+1]-1))
    {
      cumsumy.col(temp_k) = temp;
      temp = Eigen::VectorXd::Zero(ng);
      temp_k++;
    }
  }
  return cumsumy;

}

// [[Rcpp::export]]
Rcpp::List call_posindy(const Eigen::MappedSparseMatrix<double> count, const int k, const int nc)
{

  // Eigen::SparseVector<int> cck = count.col(k);
  int nnz = count.col(k).nonZeros();
  Eigen::VectorXi temp(nnz);
  Eigen::VectorXi value(nnz);
  Eigen::VectorXi ytwo(nnz);
  double mct = 0;
  int n_one = 0;
  int n_two = 0;
  int temp_k = 0;
  int nth = 0;
  Eigen::VectorXi n_onetwo(2);


  for (Eigen::MappedSparseMatrix<double>::InnerIterator it(count,k); it; ++it)
  //for (Eigen::SparseVector<int>::InnerIterator it(cck); it; ++it)
  {
    if(it.value()>0)
    {
      temp[temp_k] = it.row();
      // temp[temp_k] = it.index();
      int vt = it.value();
      value[temp_k] = vt;
      mct += vt;
      if(vt==1)
      {
        n_one++;
      }else{
        if(vt==2)
        {
          n_two++;
        }else{
          ytwo[nth++] = vt;
        }
      }
      temp_k++;
    }
  }
  // mct = mct/count.rows();
  mct = mct/nc;
  n_onetwo(0) = n_one;
  n_onetwo(1) = n_two;

  return Rcpp::List::create(Rcpp::Named("posindy") = temp.head(temp_k),
                            Rcpp::Named("Y") = value.head(temp_k),
                            Rcpp::Named("mct") = mct,
                            Rcpp::Named("n_onetwo") = n_onetwo,
                            Rcpp::Named("ytwo") = ytwo.head(nth));

}


// [[Rcpp::export]]
Rcpp::List center_m(const Eigen::Map<Eigen::MatrixXd> & X_c)
{
  Eigen::MatrixXd cm = X_c.rowwise() - X_c.colwise().mean();
  Eigen::VectorXd sds = (cm.cwiseProduct(cm)).colwise().mean();
  sds = sds.array().sqrt();
  int nc = X_c.cols();
  int nr = X_c.rows();
  for(int i=0;i<nc;i++)
  {
    if(sds(i)>0)
    {
      cm.col(i) = cm.col(i)/sds(i);
    }else{
      if(X_c(0,i)!=0)
      {
        /*
        if(X_c(0,i)!=1)
        {
          cm.col(i) = X_c.col(i)/X_c(0,i);
        }else{
          cm.col(i) = X_c.col(i);
        }
        */
        cm.col(i) = Eigen::VectorXd::Constant(nr,1);
      }else{
        sds(i) = -1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("pred") = cm,
                            Rcpp::Named("sds") = sds);

}

// [[Rcpp::export]]
Rcpp::List cv_offset(const Eigen::Map<Eigen::VectorXd> & offset_c, int input, const int nc)
{
  Eigen::VectorXd offset(nc);
  double cv = 0;
  double moffset = 1;
  if(input==1)
  {
    offset = offset_c;
    moffset = offset.mean();
  }else{
    offset = Eigen::VectorXd::Constant(nc,1);
  }
  double mexpoffset = moffset;
  
  if(moffset>0)
  {
    cv = (offset.array()-moffset).square().sum();
    cv = sqrt(cv/nc)/moffset;
  }
  
  offset = offset.array().log();
  if(input==1)
  {
    moffset = offset.mean();
  }else{
    moffset = 0;
  }
  return Rcpp::List::create(Rcpp::Named("offset") = offset,
                            Rcpp::Named("moffset") = moffset,
                            Rcpp::Named("mexpoffset") = mexpoffset,
                            Rcpp::Named("cv") = cv);

}

// [[Rcpp::export]]
double get_cv(const Eigen::Map<Eigen::VectorXd> offset_c, const Eigen::Map<Eigen::MatrixXd> X_c, 
              const Eigen::VectorXd & beta_c, const Eigen::VectorXi & cell_ind, const int ncell, const int nc)
{
  Eigen::ArrayXd extb = offset_c;
  for(int i=0; i<ncell; i++)
  {
    int ind = cell_ind(i)-1;
    extb = extb + X_c.col(ind).array()*beta_c(ind);
  }
  extb = extb.array().exp();
  double cv = 0;
  double moffset = extb.mean();
  if(moffset>0)
  {
    cv = (extb-moffset).square().sum();
    cv = cv/nc/(moffset*moffset);
  }
  return cv;
}

// [[Rcpp::export]]
Eigen::VectorXd get_cell(const Eigen::Map<Eigen::MatrixXd> X_c, const Eigen::VectorXi & fid_c, const int nb_c, const int k_c)
{
  Eigen::VectorXd iscell = Eigen::VectorXd::Zero(nb_c);
  for(int i=0; i<nb_c; i++)
  {
    for(int j=0;j<k_c;j++)
    {
      int flag = 0;
      int end = fid_c[j+1];
      int start = fid_c(j);
      double temp = X_c(start,i);
      for(int k=start;k<end;k++)
      {
        if(X_c(k,i)!=temp)
        {
          flag = 1;
          break;
        }
      }
      if(flag==1)
      {
        iscell(i) = 1;
        break;
      }
    }
  }
  return iscell;
}

// [[Rcpp::export]]
Rcpp::List ptmg_ll_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                     const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                     const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                     const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                     const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double gamma = sigma_c(1);
  double log_lambda = log(lambda);
  double log_gamma = log(gamma);

  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);

  double term1 = 0;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}


  extb = extb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;

  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::VectorXd mustar_log = mustar.log().matrix();

  term1 -= ystar.dot(mustar_log);

  term1 += k_c*alpha*log_lambda;
  term1 += nind_c*gamma*log_gamma;

  Eigen::VectorXd sum_elgcp(nind_c);
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::VectorXd gstar_phiymustar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  for(unsigned int i=0;i<nelem;i++)
  {gstar_phiymustar(posindy_c(i)) += Y_c(i);}

  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += sum_elgcp.sum();

  sum_elgcp = sum_elgcp.array() + gamma;
  gstar_phiymustar = gstar_phiymustar.array()/sum_elgcp.array();

  sum_elgcp = sum_elgcp.array().log();
  double slpey = sum_elgcp.sum();
  term1 -= gamma*slpey;
  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}
  term1 = term1*(-1);

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::MatrixXd xexb;
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  Eigen::MatrixXd dbeta_41(k_c,nb_c);

  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    xexb = X_c.block(start,0,lent,nb_c).array().colwise() * extb.segment(start,lent).array();
    xexb_f.row(i) = xexb.colwise().sum();
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb;
  }
  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
    // dbeta_42[i] = gstar_phiymustar.segment(start,lent).sum();
  }

  Eigen::ArrayXd ymumustar_dbeta_csxtb = ymumustar*(dbeta_42-cumsumxtb);
  Eigen::VectorXd db = xexb_f*ymumustar_dbeta_csxtb.matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  double dtau = 0;
  double ldm = log_lambda*k_c - mustar_log.sum();
  double adlmy = exps_s*k_c - ymustar.sum();

  dtau = dtau + alpha_pr*((cumsumxtb-dbeta_42)/mustar).sum();
  dtau = dtau + lambda_pr*ymumustar_dbeta_csxtb.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  double dtau2 = log_gamma*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return Rcpp::List::create(Rcpp::Named("fn") = term1,
                            Rcpp::Named("gr") = (-1)*der);

}


// [[Rcpp::export]]
Rcpp::List ptmg_ll_der_hes_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                             const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                             const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                             const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                             const Eigen::VectorXd & beta_c, const Eigen::VectorXd & sigma_c)
{
  double exps = exp(sigma_c(0));
  double exps_m = exps-1;
  double exps_s = sqrt(exps);
  double alpha = 1/exps_m;
  double lambda = alpha/exps_s;
  double gamma = sigma_c(1);
  double log_lambda = log(lambda);
  double log_gamma = log(gamma);

  exps_m = pow(exps_m,2);
  double alpha_pr = -exps/exps_m;
  double lambda_pr = (1-3*exps)/(2*exps_s*exps_m);
  double alpha_dpr = 2*exps*exps/(exps_m*(exps-1)) - exps/exps_m;
  double lambda_dpr = -3*exps/(2*exps_s*exps_m) + (3*exps-1)/(4*exps_s*exps_m) + (3*exps-1)*exps_s/(exps_m*(exps-1));

  double term1 = 0;
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nb_c+2,nb_c+2);
  double hes_sigma = alpha_dpr*log_lambda + 2*alpha_pr*lambda_pr/lambda - alpha*pow(lambda_pr/lambda,2) + alpha/lambda*lambda_dpr;
  hes_sigma = hes_sigma*k_c;

  Eigen::VectorXd extb = offset_c + X_c*beta_c;
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {term1 += extb(posindy_c(i))*Y_c(i);}

  extb = extb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
    cumsumxtb(i)= extb.segment(fid_c(i),len(i)).sum();
  }

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;

  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::VectorXd mustar_log = mustar.log().matrix();

  term1 -= ystar.dot(mustar_log);

  term1 += k_c*alpha*log_lambda;
  term1 += nind_c*gamma*log_gamma;

  Eigen::VectorXd sum_elgcp(nind_c);
  Eigen::ArrayXd ymustar = ystar.array()/mustar;
  Eigen::ArrayXd ymumustar = ymustar/mustar;
  Eigen::VectorXd gstar_phiymustar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  for(unsigned int i=0;i<nelem;i++)
  {gstar_phiymustar(posindy_c(i)) += Y_c(i);}

  for(int i=0;i<k_c;i++)
  {
    sum_elgcp.segment(fid_c(i),len(i)) = ymustar(i)*extb.segment(fid_c(i),len(i));
  }

  term1 += sum_elgcp.sum();

  sum_elgcp = sum_elgcp.array() + gamma;
  Eigen::ArrayXd tempa = 1/sum_elgcp.array();
  hessian(nb_c+1,nb_c+1) = nind_c/gamma - 2*tempa.sum();
  gstar_phiymustar = gstar_phiymustar.array()*tempa;

  sum_elgcp = sum_elgcp.array().log();
  double slpey = sum_elgcp.sum();
  term1 -= gamma*slpey;
  for(unsigned int i=0;i<nelem;i++)
  {term1 -= Y_c(i)*sum_elgcp(posindy_c(i));}
  term1 = term1*(-1);

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c+2);

  Eigen::MatrixXd xexb = X_c.array().colwise()*extb.array();
  Eigen::MatrixXd xexb_f(k_c,nb_c);
  Eigen::MatrixXd dbeta_41(k_c,nb_c);

  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // xexb = X_c.block(start,0,lent,nb_c).array().colwise() * extb.segment(start,lent).array();
    xexb_f.row(i) = xexb.block(start,0,lent,nb_c).colwise().sum();
    dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
  }
  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd dbeta_42(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start,lent).dot(extb.segment(start,lent));
    // dbeta_42[i] = gstar_phiymustar.segment(start,lent).sum();
  }

  dbeta_42 = dbeta_42-cumsumxtb;

  Eigen::ArrayXd ymumustar_dbeta_csxtb = ymumustar*dbeta_42;
  Eigen::VectorXd db = xexb_f*ymumustar_dbeta_csxtb.matrix() - dbeta_41*ymustar.matrix();
  for(unsigned int i=0;i<nelem;i++)
  {db += X_c.row(posindy_c(i)).transpose()*Y_c(i);}
  double dtau = 0;
  double ldm = log_lambda*k_c - mustar_log.sum();
  double adlmy = exps_s*k_c - ymustar.sum();

  Eigen::ArrayXd imustar = 1/mustar;
  Eigen::VectorXd hes_b_l = dbeta_42*imustar;
  double dbim = hes_b_l.sum();
  dtau = dtau - alpha_pr*dbim;
  dtau = dtau + lambda_pr*ymumustar_dbeta_csxtb.sum();
  dtau += alpha_pr*ldm + lambda_pr*adlmy;

  // not necessary
  hes_sigma -= dbim*alpha_dpr - 2*alpha_pr*lambda_pr*hes_b_l.dot(imustar.matrix()) - lambda_dpr*ymumustar_dbeta_csxtb.sum() + 2*lambda_pr*lambda_pr*(ymumustar_dbeta_csxtb*imustar).sum();

  double dtau2 = log_gamma*nind_c + nind_c - slpey - gstar_phiymustar.sum();

  // Eigen::VectorXd hes_sigma_beta = -alpha_pr*(xexb_f*imustar.matrix()) + lambda_pr*(xexb_f*ymumustar.matrix());
  Eigen::VectorXd hes_sigma_beta = -alpha_pr*(dbeta_41*imustar.matrix()) + lambda_pr*(dbeta_41*ymumustar.matrix());
  hes_sigma_beta -= -alpha_pr*(xexb_f*hes_b_l.cwiseProduct(imustar.matrix())) + 2*lambda_pr*(xexb_f*hes_b_l.cwiseProduct(ymumustar.matrix()));

  hes_sigma = hes_sigma - (2*alpha_pr*lambda_pr*imustar.sum() + alpha_dpr*mustar_log.sum() + lambda_dpr*ymustar.sum() - lambda_pr*lambda_pr*ymumustar.sum());
  Eigen::VectorXd apmmlpymm = alpha_pr*imustar - lambda_pr*ymumustar;


  hes_b_l = gstar_phiymustar.array()*tempa;
  hessian(nb_c+1,nb_c+1) = hessian(nb_c+1,nb_c+1) + hes_b_l.sum();

  Eigen::VectorXd hes_b_l_f(k_c);
  Eigen::MatrixXd hes_b_l_xexb = xexb.array().colwise()*(hes_b_l.array()-tempa);
  Eigen::MatrixXd hes_b_l_xexb_f(k_c,nb_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start,0,lent,nb_c).colwise().sum();
    // hes_b_l_f(i) = hes_b_l.segment(start,lent).sum();
  }
  Eigen::VectorXd hes_gamma_beta = hes_b_l_xexb_f.transpose()*ymustar.matrix();

  hes_b_l = hes_b_l.array()*extb.array();
  Eigen::VectorXd hes_b2_l = hes_b_l - extb.cwiseProduct(tempa.matrix());
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_f(i) = hes_b2_l.segment(start,lent).sum();
  }
  hes_gamma_beta = hes_gamma_beta - xexb_f*hes_b_l_f.cwiseProduct(ymumustar.matrix());
  hessian.row(nb_c+1).head(nb_c) = hes_gamma_beta;
  hessian.col(nb_c+1).head(nb_c) = hes_gamma_beta;
  hessian(nb_c+1,nb_c) = apmmlpymm.dot(hes_b_l_f);
  hessian(nb_c,nb_c+1) = hessian(nb_c+1,nb_c);

  hes_b2_l = hes_b_l.array()*extb.array();

  Eigen::VectorXd tempb;
  Eigen::VectorXd tempc;
  hes_b_l_xexb = xexb.array().colwise()*hes_b_l.array();
  Eigen::VectorXd hes_b2_l_f(k_c);
  for(int i=0;i<k_c;i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start,0,lent,nb_c).colwise().sum();
    hes_b2_l_f(i) = hes_b2_l.segment(start,lent).sum();
  }
  // hes_b_l_xexb_f.transposeInPlace();

  hes_sigma_beta += hes_b_l_xexb_f.transpose()*(apmmlpymm.cwiseProduct(ymustar.matrix())) - xexb_f*(apmmlpymm.cwiseProduct(hes_b2_l_f.cwiseProduct(ymumustar.matrix())));
  hessian.row(nb_c).head(nb_c) = hes_sigma_beta;
  hessian.col(nb_c).head(nb_c) = hes_sigma_beta;

  apmmlpymm = apmmlpymm.cwiseProduct(apmmlpymm);
  hes_sigma += hes_b2_l_f.dot(apmmlpymm);
  hessian(nb_c,nb_c) = hes_sigma;

  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();
  apmmlpymm = ymustar.cwiseProduct(ymustar);
  Eigen::VectorXd ymuymuimu = apmmlpymm.cwiseProduct(imustar.matrix());
  hes_b2_l_f = hes_b2_l_f.array()*ymuymuimu.array()*imustar;

  for(int i=0;i<nb_c;i++)
  {
    for(int j=i;j<nb_c;j++)
    {
      tempb = xexb.col(i).cwiseProduct(X_c.col(j));
      tempc = tempb.cwiseProduct(hes_b_l);

      for(int k=0;k<k_c;k++)
      {
        int start = fid_c[k];
        int lent = len(k);

        hessian(i,j) -= tempb.segment(start,lent).dot(gstar_phiymustar.segment(start,lent))*ymustar(k);
        hessian(i,j) += tempc.segment(start,lent).sum()*apmmlpymm(k);
        //not necessary
        hessian(i,j) += ymumustar_dbeta_csxtb(k)*tempb.segment(start,lent).sum();
      }
      hes_b_l_f = dbeta_41.col(i).cwiseProduct(xexb_f.col(j)) + dbeta_41.col(j).cwiseProduct(xexb_f.col(i));
      hessian(i,j) += hes_b_l_f.dot(ymumustar.matrix());
      hessian(i,j) -= ymuymuimu.dot(hes_b_l_xexb_f.col(i).cwiseProduct(xexb_f.col(j))+hes_b_l_xexb_f.col(j).cwiseProduct(xexb_f.col(i)));
      tempb = xexb_f.col(i).cwiseProduct(xexb_f.col(j));
      hessian(i,j) -= tempb.dot(ymumustar.matrix());
      hessian(i,j) += hes_b2_l_f.dot(tempb);
      //not necessary
      hessian(i,j) += (-2)*tempb.dot((ymumustar_dbeta_csxtb*imustar).matrix());

      if(i!=j)
      {hessian(j,i) = hessian(i,j);}
    }
  }

  der.segment(0,nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c+1] = dtau2;

  return Rcpp::List::create(Rcpp::Named("fn") = term1,
                            Rcpp::Named("gr") = (-1)*der,
                            Rcpp::Named("hes") = (-1)*hessian);

}


// [[Rcpp::export]]
Rcpp::List opt_pml_hl(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                             const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                             const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                             const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                             const Eigen::VectorXd & beta_c,
                             const Eigen::VectorXd & sigma_c, const int reml, const double eps, const int ord)
{
  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double loglik = 0;

  Eigen::VectorXd logw = Eigen::VectorXd::Constant(k_c,1,0);
  // Eigen::VectorXd logw = logw_c;
  Eigen::VectorXd beta = beta_c;

  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}

  Eigen::VectorXd y_full = Eigen::VectorXd::Zero(nind_c);
  for(unsigned int i=0; i<nelem; i++)
  {
    int idx = posindy_c(i);
    if(idx >= 0 && idx < nind_c)
    {
      y_full(idx) = Y_c(i);
    }
  }

  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }

  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  Eigen::VectorXd extb = offset_c + X_c*beta;

  Eigen::VectorXd w = logw.array().exp();

  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}

  loglik += logw.dot(cumsumy_c);

  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();

  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}

  loglik += alpha*logw.sum() - lambda*w.sum();
  
  //double eps = 1e-6;
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);
  int stepd = 0;
  int maxstd = 10;
  int maxstep = 50;
  double convd = 0.01;

  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||((likdif>eps)&&(step<maxstep)))
  {
    step++;

    // double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);
    Eigen::VectorXd damp = Eigen::VectorXd::Constant(nb_c,1);

    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
    Eigen::VectorXd dw(k_c);
    for(int i=0;i<k_c;i++)
    {
      dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    dw = cumsumy_c - dw - lambda*w;
    dw = dw.array() + alpha;

    gstar_extb_phi = gstar_extb_phi.array()/(extb.array()+gamma);

    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;

    Eigen::MatrixXd xgsetbp = X_c.array().colwise()*gstar_extb_phi.array();

    for(int i=0;i<k_c;i++)
    {
      vwb.row(i) = xgsetbp.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }
    vwb = gamma*vwb;
    for(int i=0;i<nb_c;i++)
    {
      for(int j=i;j<nb_c;j++)
      {
        vb(i,j) = X_c.col(i).dot(xgsetbp.col(j));
        if(i!=j)
        {
          vb(j,i) = vb(i,j);
        }
      }
    }
    // vb = X_c.transpose()*vb;
    vb = gamma*vb;
    //Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    //vb2 = vb - vwb.transpose()*temp;
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array().sqrt();
    vb2 = vb - temp.transpose()*temp;

    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    //beta = beta + stepbeta;
    Eigen::VectorXd new_b = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    //logw = logw + steplogw;
    Eigen::VectorXd new_w = logw + steplogw;

    loglikp = loglik;
    loglik = 0;

    extb = offset_c + X_c*new_b;
    w = new_w.array().exp();

    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += new_w.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
    }
    extb = extb.array().exp();

    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik += alpha*new_w.sum() - lambda*w.sum();

    likdif = loglik - loglikp;
    stepd = 0;
    double minstep = 40;
    
    while((likdif<0)||(std::isinf(loglik)))
    {
      stepd++;
      minstep = minstep/2;
      
      if(stepd>maxstd)
      {
        likdif = 0;
        loglik = loglikp;
        double mabsdb = db.cwiseAbs().maxCoeff();
        double mabsdw = dw.cwiseAbs().maxCoeff();
        if((mabsdb>convd) || (mabsdw>convd))
        {stepd++;}
        break;
      }
      
      for(int i=0;i<nb_c;i++)
      {
        if((stepbeta(i)<40)&&(stepbeta(i)>-40))
          {
            damp(i) = damp(i)/2;
            new_b(i) = beta(i) + stepbeta(i)*damp(i);
          }else{
              if(stepbeta(i)>0)
              {
                new_b(i) = beta(i) + minstep;
              }else{
                new_b(i) = beta(i) - minstep;
              }
          }
      }
      
      for(int i=0;i<k_c;i++)
      {
        if((steplogw(i)<40)&&(steplogw(i)>-40))
          {
            damp_w(i)=damp_w(i)/2;
            new_w(i) = logw(i) + steplogw(i)*damp_w(i);
          }else{
            if(steplogw(i)>0)
            {
              new_w(i) = logw(i) + minstep;
            }else{
              new_w(i) = logw(i) - minstep;
            }
          }
      }

      loglik = 0;

      extb = offset_c + X_c*new_b;
      w = new_w.array().exp();

      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += new_w.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
      }
      extb = extb.array().exp();

      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik += alpha*new_w.sum() - lambda*w.sum();

      likdif = loglik - loglikp;

    }
    beta = new_b;
    logw = new_w;
  }
  
  double logdet;
  if(reml==1)
  {
    logdet = vw.array().abs().log().sum();
    double absdet = vb2.determinant();
    if(absdet<0)
      absdet *= -1;
    logdet += log(absdet);
  }else{
    logdet = vw.array().abs().log().sum();
  }
  
  double sec_ord = 0;
  Eigen::VectorXd sec_ord_vec = Eigen::VectorXd::Zero(k_c);
  if(ord>1)
  {
    Eigen::VectorXd third_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd four_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd temp_der = Eigen::VectorXd::Zero(k_c);
    
    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::ArrayXd extbg = extb.array()+gamma;
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;
    Eigen::ArrayXd vws = vw.array()*vw.array();
    
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    
    gstar = gstar_extb_phi.array()*(gamma-extb.array());
    for(int i=0;i<k_c;i++)
    {
      third_der(i) = gstar.segment(fid_c(i),len(i)).sum();
    }
    third_der = gamma*third_der;
    third_der += lambda*w;
    temp_der = third_der.array()*third_der.array()/(vws*vw.array());
    sec_ord_vec += (5.0/24.0)*temp_der.matrix();
    sec_ord += 5*temp_der.sum()/24;

    if(ord>2)
    {

      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      Eigen::ArrayXd extbp = extb.array()*extb.array();
      gstar = gstar_extb_phi.array()*(gamma*gamma+extbp-(4*gamma)*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()/vws;
      sec_ord_vec -= (1.0/8.0)*temp_der.matrix();
      sec_ord -= temp_der.sum()/8;

      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      gstar = gstar_extb_phi.array()*(gamma*gamma*gamma-(11*gamma*gamma)*extb.array()+(11*gamma)*extbp-extbp*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()*third_der.array()/(vws*vws);
      sec_ord_vec += (7.0/48.0)*temp_der.matrix();
      sec_ord += 7*temp_der.sum()/48;
    }

  }

  double phi = 1.0 / gamma;
  double inv_phi = gamma;
  double phi_sq = phi * phi;
  double kappa_s = 1.0 / (std::exp(sigma_c(0)) - 1.0);
  double digamma_kappa = R::digamma(kappa_s);

  Eigen::MatrixXd h_beta = Eigen::MatrixXd::Zero(nb_c, k_c);
  Eigen::MatrixXd h_eta_beta = Eigen::MatrixXd::Zero(nb_c, k_c);
  Eigen::MatrixXd h_etaeta_beta = Eigen::MatrixXd::Zero(nb_c, k_c);
  Rcpp::List h_beta_beta_list(k_c);
  Eigen::MatrixXd h_beta_phi = Eigen::MatrixXd::Zero(nb_c, k_c);
  Eigen::VectorXd h_phi = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_eta_phi = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_etaeta_phi = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_phi_phi = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_tau = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_eta_tau = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_etaeta_tau = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_tau_tau = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_etaeta_vec = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_etaetaeta_vec = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd h_etaetaetaeta_vec = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd log_H_etaeta = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd loglik_subject = Eigen::VectorXd::Zero(k_c);
  Eigen::VectorXd log_vw_vec = Eigen::VectorXd::Zero(k_c);
  Rcpp::List h_etaeta_beta_beta_list(k_c);
  Rcpp::List g_beta_beta_list(k_c);
  Rcpp::List h_etaeta_beta_phi_list(k_c);
  Rcpp::List f_beta_phi_list(k_c);
  Rcpp::List g_beta_phi_list(k_c);
  Rcpp::NumericVector h_etaeta_phi_phi_vec(k_c);
  Rcpp::NumericVector g_phi_phi_vec(k_c);
  Rcpp::NumericVector b_phi_vec(k_c);
  Eigen::MatrixXd b_beta = Eigen::MatrixXd::Zero(nb_c, k_c);

  double dkappa_dtau = -exps / ((exps - 1.0) * (exps - 1.0));
  double trigamma_kappa = R::trigamma(kappa_s);

  for(int subj = 0; subj < k_c; ++subj)
  {
    int start = fid_c[subj];
    int end = fid_c[subj + 1];
    int len = end - start;

    Eigen::MatrixXd Xi = X_c.block(start, 0, len, nb_c);
    Eigen::VectorXd offset_i = offset_c.segment(start, len);
    Eigen::VectorXd y_i = y_full.segment(start, len);
    Eigen::VectorXd base_lin = offset_i + Xi * beta;

    double eta_i = logw(subj);
    Eigen::VectorXd lin = offset_i + Xi * beta;
    lin.array() += eta_i;
    Eigen::ArrayXd mu = lin.array().exp();
    Eigen::ArrayXd phi_mu = phi * mu;
    Eigen::ArrayXd denom = 1.0 + phi_mu;
    Eigen::ArrayXd y_plus = y_i.array() + inv_phi;

    Eigen::VectorXd h_beta_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::VectorXd h_eta_beta_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::VectorXd h_etaeta_beta_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::MatrixXd h_etaeta_beta_beta_i = Eigen::MatrixXd::Zero(nb_c, nb_c);
    Eigen::MatrixXd h_beta_beta_i = Eigen::MatrixXd::Zero(nb_c, nb_c);
    Eigen::MatrixXd g_beta_beta_i = Eigen::MatrixXd::Zero(nb_c, nb_c);
    Eigen::VectorXd h_etaeta_beta_phi_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::VectorXd f_beta_phi_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::VectorXd g_beta_phi_i = Eigen::VectorXd::Zero(nb_c);
    Eigen::VectorXd b_beta_i = Eigen::VectorXd::Zero(nb_c);
    double h_phi_i = 0.0;
    double h_eta_phi_i = 0.0;
    double h_etaeta_phi_i = 0.0;
    double h_etaeta_phi_phi_i = 0.0;
    double h_phi_phi_i = 0.0;
    double g_phi_phi_i = 0.0;
    double b_phi_i = 0.0;
    double sum_t = 0.0;
    double sum_u = 0.0;
    double sum_v = 0.0;
    Eigen::ArrayXd log_mu_plus_gamma(len);

    for(int r = 0; r < len; ++r)
    {
      double y = y_i(r);
      double mu_r = mu(r);
      double phi_mu_r = phi_mu(r);
      double denom_r = denom(r);
      double y_plus_r = y_plus(r);
      double log_denom_r = std::log(denom_r);

      double s = y - y_plus_r * (phi_mu_r / denom_r);
      double t = - y_plus_r * (phi_mu_r / (denom_r * denom_r));
      double u = - y_plus_r * (phi_mu_r * (1.0 - phi_mu_r) / (denom_r * denom_r * denom_r));
      double v = - y_plus_r * (phi_mu_r) * (1.0 - 4.0 * phi_mu_r + phi_mu_r * phi_mu_r) /
        (denom_r * denom_r * denom_r * denom_r);

      double A = y_plus_r;
      double Aprime = -inv_phi * inv_phi;
      double Aprime2 = 2.0 * inv_phi * inv_phi * inv_phi;
      double B = phi * mu_r;
      double Bprime = mu_r;
      double C = 1.0 - B;
      double Cprime = -Bprime;
      double D = 1.0 + B;
      double Dprime = Bprime;

      double f_t = -A * B;
      double g_t = 1.0 / (D * D);
      double f_t_prime = -(Aprime * B + A * Bprime);
      double g_t_prime = -2.0 * Dprime / (D * D * D);
      double f_t_prime2 = -(Aprime2 * B + 2.0 * Aprime * Bprime);
      double g_t_prime2 = 6.0 * Dprime * Dprime / (D * D * D * D);

      double t_phi = f_t_prime * g_t + f_t * g_t_prime;
      double t_phi_phi = f_t_prime2 * g_t + 2.0 * f_t_prime * g_t_prime + f_t * g_t_prime2;

      double f_u = -A * B * C;
      double g_u = 1.0 / (D * D * D);
      double f_u_prime = -(Aprime * B * C + A * Bprime * C + A * B * Cprime);
      double g_u_prime = -3.0 * Dprime / (D * D * D * D);
      double u_phi = f_u_prime * g_u + f_u * g_u_prime;

      double digamma_term = R::digamma(y + inv_phi) - R::digamma(inv_phi);
      double trigamma_diff = R::trigamma(y + inv_phi) - R::trigamma(inv_phi);
      double ell_phi_phi = -y / (phi * phi)
        + (2.0 / (phi * phi * phi)) * digamma_term
        + (1.0 / (phi * phi * phi * phi)) * trigamma_diff
        - (2.0 / (phi * phi * phi)) * std::log(denom_r)
        + (2.0 / (phi * phi)) * (mu_r / denom_r)
        + (y + inv_phi) * (mu_r * mu_r) / (denom_r * denom_r);

      h_beta_i += Xi.row(r).transpose() * s;
      h_eta_beta_i += Xi.row(r).transpose() * t;
      h_etaeta_beta_i += Xi.row(r).transpose() * u;
      h_beta_beta_i += Xi.row(r).transpose() * Xi.row(r) * t;
      h_etaeta_beta_beta_i += Xi.row(r).transpose() * Xi.row(r) * u;
      g_beta_beta_i += Xi.row(r).transpose() * Xi.row(r) * v;
      b_beta_i += Xi.row(r).transpose() * v;
      h_etaeta_beta_phi_i += Xi.row(r).transpose() * u_phi;
      f_beta_phi_i += Xi.row(r).transpose() * t_phi;
      g_beta_phi_i += Xi.row(r).transpose() * u_phi;

      double ell_phi = (y / phi) - (digamma_term / (phi_sq)) + (log_denom_r / (phi_sq)) - y_plus_r * (mu_r / denom_r);
      double s_phi = (phi_mu_r / denom_r) / (phi_sq) - y_plus_r * (mu_r / (denom_r * denom_r));

      h_phi_i += ell_phi;
      h_eta_phi_i += s_phi;
      h_etaeta_phi_i += t_phi;
      h_beta_phi.col(subj) += Xi.row(r).transpose() * s_phi;
      h_phi_phi_i += ell_phi_phi;
      h_etaeta_phi_phi_i += t_phi_phi;
      g_phi_phi_i += t_phi_phi;
      b_phi_i += u_phi;

      sum_t += t;
      sum_u += u;
      sum_v += v;

      log_mu_plus_gamma(r) = std::log(mu_r + gamma);
    }

    double exp_eta = std::exp(eta_i);
    double h_etaeta_prior = -kappa_s * exp_eta;
    double h_etaetaeta_prior = -kappa_s * exp_eta;
    double h_eta_tau_i_val = -kappa_s * (1.0 - exp_eta);
    double h_etaeta_tau_i_val = kappa_s * exp_eta;
    double h_tau_i_val = -kappa_s * (eta_i - exp_eta + std::log(kappa_s) + 1.0 - digamma_kappa);

    double h_etaeta_i = sum_t + h_etaeta_prior;
    double h_etaetaeta_i = sum_u + h_etaetaeta_prior;
    double h_etaetaetaeta_i = sum_v + h_etaetaeta_prior;

    double H_val = -h_etaeta_i;
    if (H_val <= 0.0) {
      H_val = 1e-12;
    }

    double logw_i = eta_i;
    double w_i = std::exp(logw_i);
    double vw_i = vw(subj);
    if (vw_i <= 0.0) {
      vw_i = 1e-12;
    }

    double loglik_i = 0.0;
    for (int r = 0; r < len; ++r) {
      double y = y_i(r);
      if (y > 0) {
        loglik_i += base_lin(r) * y;
      }
    }
    loglik_i += logw_i * cumsumy_c(subj);
    loglik_i -= gamma * log_mu_plus_gamma.sum();
    loglik_i -= (y_i.array() * log_mu_plus_gamma).sum();
    loglik_i += alpha * logw_i - lambda * w_i;

    h_beta.col(subj) = h_beta_i;
    h_eta_beta.col(subj) = h_eta_beta_i;
    h_etaeta_beta.col(subj) = h_etaeta_beta_i;
    h_phi(subj) = h_phi_i;
    h_eta_phi(subj) = h_eta_phi_i;
    h_etaeta_phi(subj) = h_etaeta_phi_i;
    h_tau(subj) = h_tau_i_val;
    h_eta_tau(subj) = h_eta_tau_i_val;
    h_etaeta_tau(subj) = h_etaeta_tau_i_val;
    h_etaeta_vec(subj) = h_etaeta_i;
    h_etaetaeta_vec(subj) = h_etaetaeta_i;
    h_etaetaetaeta_vec(subj) = h_etaetaetaeta_i;
    h_beta_beta_list[subj] = Rcpp::wrap(h_beta_beta_i);
    h_etaeta_beta_beta_list[subj] = Rcpp::wrap(h_etaeta_beta_beta_i);
    g_beta_beta_list[subj] = Rcpp::wrap(g_beta_beta_i);
    h_etaeta_beta_phi_list[subj] = Rcpp::wrap(h_etaeta_beta_phi_i);
    f_beta_phi_list[subj] = Rcpp::wrap(f_beta_phi_i);
    g_beta_phi_list[subj] = Rcpp::wrap(g_beta_phi_i);
    b_beta.col(subj) = b_beta_i;
    h_etaeta_phi_phi_vec[subj] = h_etaeta_phi_phi_i;
    h_phi_phi(subj) = h_phi_phi_i;
    g_phi_phi_vec[subj] = g_phi_phi_i;
    b_phi_vec[subj] = b_phi_i;
    log_H_etaeta(subj) = std::log(H_val);
    loglik_subject(subj) = loglik_i;
    log_vw_vec(subj) = std::log(vw_i);

    double term_bracket = eta_i - std::exp(eta_i) + std::log(kappa_s) + 1.0 - digamma_kappa;
    h_tau_tau(subj) = dkappa_dtau * (-term_bracket - 1.0 + kappa_s * trigamma_kappa);
  }

  Rcpp::List per_subject = Rcpp::List::create(
    Rcpp::Named("h_beta") = h_beta,
    Rcpp::Named("h_eta_beta") = h_eta_beta,
    Rcpp::Named("h_etaeta_beta") = h_etaeta_beta,
    Rcpp::Named("h_beta_beta") = h_beta_beta_list,
    Rcpp::Named("h_beta_phi") = h_beta_phi,
    Rcpp::Named("h_phi") = h_phi,
    Rcpp::Named("h_eta_phi") = h_eta_phi,
    Rcpp::Named("h_etaeta_phi") = h_etaeta_phi,
    Rcpp::Named("h_phi_phi") = h_phi_phi,
    Rcpp::Named("h_tau") = h_tau,
    Rcpp::Named("h_eta_tau") = h_eta_tau,
    Rcpp::Named("h_etaeta_tau") = h_etaeta_tau,
    Rcpp::Named("h_tau_tau") = h_tau_tau,
    Rcpp::Named("h_etaeta") = h_etaeta_vec,
    Rcpp::Named("h_etaetaeta") = h_etaetaeta_vec,
    Rcpp::Named("h_etaetaetaeta") = h_etaetaetaeta_vec,
    Rcpp::Named("h_etaeta_beta_beta") = h_etaeta_beta_beta_list,
    Rcpp::Named("g_beta_beta") = g_beta_beta_list,
    Rcpp::Named("h_etaeta_beta_phi") = h_etaeta_beta_phi_list,
    Rcpp::Named("f_beta_phi") = f_beta_phi_list,
    Rcpp::Named("g_beta_phi") = g_beta_phi_list,
    Rcpp::Named("h_etaeta_phi_phi") = h_etaeta_phi_phi_vec,
    Rcpp::Named("g_phi_phi") = g_phi_phi_vec,
    Rcpp::Named("b_phi") = b_phi_vec,
    Rcpp::Named("g_tau_tau") = (-1.0) * h_etaeta_tau,
    Rcpp::Named("b_beta") = b_beta,
    Rcpp::Named("log_H_etaeta") = log_H_etaeta,
    Rcpp::Named("loglik") = loglik_subject,
    Rcpp::Named("sec_ord") = sec_ord_vec,
    Rcpp::Named("log_vw") = log_vw_vec
  );

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            Rcpp::Named("var") = vb2,
                            Rcpp::Named("h_etaeta") = vw,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("damp") = stepd,
                            Rcpp::Named("second") = sec_ord,
                            Rcpp::Named("per_subject_stats") = per_subject);

}

// [[Rcpp::export]]
Rcpp::List opt_pml(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                             const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                             const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                             const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                             const Eigen::VectorXd & beta_c,
                             const Eigen::VectorXd & sigma_c, const int reml, const double eps, const int ord)
{
  double exps = exp(sigma_c(0));
  double alpha = 1/(exps-1);
  double lambda = 1/(sqrt(exps)*(exps-1));
  double gamma = sigma_c(1);

  double loglik = 0;

  Eigen::VectorXd logw = Eigen::VectorXd::Constant(k_c,1,0);
  Eigen::VectorXd beta = beta_c;

  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}

  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }

  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  Eigen::VectorXd extb = offset_c + X_c*beta;

  Eigen::VectorXd w = logw.array().exp();

  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}

  loglik += logw.dot(cumsumy_c);

  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();

  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}

  loglik += alpha*logw.sum() - lambda*w.sum();
  
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);
  int stepd = 0;
  int maxstd = 10;
  int maxstep = 50;
  double convd = 0.01;

  while((step==0)||((likdif>eps)&&(step<maxstep)))
  {
    step++;

    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);
    Eigen::VectorXd damp = Eigen::VectorXd::Constant(nb_c,1);

    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
    Eigen::VectorXd dw(k_c);
    for(int i=0;i<k_c;i++)
    {
      dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    dw = cumsumy_c - dw - lambda*w;
    dw = dw.array() + alpha;

    gstar_extb_phi = gstar_extb_phi.array()/(extb.array()+gamma);

    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;

    Eigen::MatrixXd xgsetbp = X_c.array().colwise()*gstar_extb_phi.array();

    for(int i=0;i<k_c;i++)
    {
      vwb.row(i) = xgsetbp.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }
    vwb = gamma*vwb;
    for(int i=0;i<nb_c;i++)
    {
      for(int j=i;j<nb_c;j++)
      {
        vb(i,j) = X_c.col(i).dot(xgsetbp.col(j));
        if(i!=j)
        {
          vb(j,i) = vb(i,j);
        }
      }
    }
    vb = gamma*vb;
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array().sqrt();
    vb2 = vb - temp.transpose()*temp;

    Eigen::VectorXd dwvw = dw.array()/vw.array();
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd new_b = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    Eigen::VectorXd new_w = logw + steplogw;

    loglikp = loglik;
    loglik = 0;

    extb = offset_c + X_c*new_b;
    w = new_w.array().exp();

    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += new_w.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
    }
    extb = extb.array().exp();

    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik += alpha*new_w.sum() - lambda*w.sum();

    likdif = loglik - loglikp;

    if((likdif<0)||(std::isinf(loglik)))
    {
      stepd++;
      double minstep = 40;
      
      while((likdif<0)||(std::isinf(loglik)))
      {
        minstep = minstep/2;
        
        if(stepd>maxstd)
        {
          likdif = 0;
          loglik = loglikp;
          double mabsdb = db.cwiseAbs().maxCoeff();
          double mabsdw = dw.cwiseAbs().maxCoeff();
          if((mabsdb>convd) || (mabsdw>convd))
          {stepd++;}
          break;
        }
        
        for(int i=0;i<nb_c;i++)
        {
          if((stepbeta(i)<40)&&(stepbeta(i)>-40))
          {
            damp(i) = damp(i)/2;
            new_b(i) = beta(i) + stepbeta(i)*damp(i);
          }else{
            if(stepbeta(i)>0)
            {
              new_b(i) = beta(i) + minstep;
            }else{
              new_b(i) = beta(i) - minstep;
            }
          }
        }
        
        for(int i=0;i<k_c;i++)
        {
          if((steplogw(i)<40)&&(steplogw(i)>-40))
          {
            damp_w(i)=damp_w(i)/2;
            new_w(i) = logw(i) + steplogw(i)*damp_w(i);
          }else{
            if(steplogw(i)>0)
            {
              new_w(i) = logw(i) + minstep;
            }else{
              new_w(i) = logw(i) - minstep;
            }
          }
        }

        loglik = 0;

        extb = offset_c + X_c*new_b;
        w = new_w.array().exp();

        for(unsigned int i=0;i<nelem;i++)
        {loglik += extb(posindy_c(i))*Y_c(i);}
        loglik += new_w.dot(cumsumy_c);
        for(int i=0;i<k_c;i++)
        {
          extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
        }
        extb = extb.array().exp();

        extbphil = (extb.array()+gamma).log();
        loglik -= gamma*extbphil.sum();
        for(unsigned int i=0;i<nelem;i++)
        {loglik -= Y_c(i)*extbphil(posindy_c(i));}
        loglik += alpha*new_w.sum() - lambda*w.sum();

        likdif = loglik - loglikp;

        stepd++;
      }
    }
    beta = new_b;
    logw = new_w;
  }
  
  double logdet;
  if(reml==1)
  {
    logdet = vw.array().abs().log().sum();
    double absdet = vb2.determinant();
    if(absdet<0)
      absdet *= -1;
    logdet += log(absdet);
  }else{
    logdet = vw.array().abs().log().sum();
  }
  
  double sec_ord = 0;
  if(ord>1)
  {
    Eigen::VectorXd third_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd four_der = Eigen::VectorXd::Zero(k_c);
    Eigen::VectorXd temp_der = Eigen::VectorXd::Zero(k_c);
    
    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::ArrayXd extbg = extb.array()+gamma;
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw += lambda*w;
    Eigen::ArrayXd vws = vw.array()*vw.array();
    
    gstar_extb_phi = gstar_extb_phi.array()/extbg;
    
    gstar = gstar_extb_phi.array()*(gamma-extb.array());
    for(int i=0;i<k_c;i++)
    {
      third_der(i) = gstar.segment(fid_c(i),len(i)).sum();
    }
    third_der = gamma*third_der;
    third_der += lambda*w;
    temp_der = third_der.array()*third_der.array()/(vws*vw.array());
    sec_ord += 5*temp_der.sum()/24;
    
    if(ord>2)
    {
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      Eigen::ArrayXd extbp = extb.array()*extb.array();
      gstar = gstar_extb_phi.array()*(gamma*gamma+extbp-(4*gamma)*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()/vws;
      sec_ord -= temp_der.sum()/8;
      
      gstar_extb_phi = gstar_extb_phi.array()/extbg;
      gstar = gstar_extb_phi.array()*(gamma*gamma*gamma-(11*gamma*gamma)*extb.array()+(11*gamma)*extbp-extbp*extb.array());
      for(int i=0;i<k_c;i++)
      {
        four_der(i) = gstar.segment(fid_c(i),len(i)).sum();
      }
      four_der = gamma*four_der;
      four_der += lambda*w;
      temp_der = four_der.array()*third_der.array()/(vws*vws);
      sec_ord += 7*temp_der.sum()/48;
    }
     
  }

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            Rcpp::Named("var") = vb2,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("damp") = stepd,
                            Rcpp::Named("second") = sec_ord);

}

// [[Rcpp::export]]
Rcpp::List opt_pml_constrained(
    const Eigen::Map<Eigen::MatrixXd> & X_c,
    const Eigen::Map<Eigen::VectorXd> & offset_c,
    const Eigen::VectorXd & Y_c,
    const Eigen::VectorXi & fid_c,
    const Eigen::VectorXd & cumsumy_c,
    const Eigen::VectorXi & posind_c,
    const Eigen::VectorXi & posindy_c,
    const Eigen::VectorXd & beta_anchor,
    const Eigen::Map<Eigen::MatrixXd> & K_c,
    const Eigen::VectorXd & gamma_init,
    const Eigen::VectorXd & sigma_c,
    const int reml,
    const double eps,
    const int ord) {

  const int free_dim = K_c.cols();
  const int subject_count = fid_c.size() - 1;
  if (free_dim <= 0) {
    return Rcpp::List::create(Rcpp::Named("beta") = beta_anchor,
                              Rcpp::Named("var") = Eigen::MatrixXd::Zero(beta_anchor.size(), beta_anchor.size()),
                              Rcpp::Named("loglik") = 0.0,
                              Rcpp::Named("loglikp") = 0.0,
                              Rcpp::Named("h_etaeta") = Eigen::VectorXd::Zero(fid_c.size() - 1),
                              Rcpp::Named("logdet") = 0.0,
                              Rcpp::Named("iter") = 0,
                              Rcpp::Named("damp") = 0,
                              Rcpp::Named("second") = 0.0,
                              Rcpp::Named("logw") = Eigen::VectorXd::Zero(fid_c.size() - 1));
  }

  Eigen::MatrixXd X_gamma = X_c * K_c;
  Eigen::VectorXd offset_gamma = offset_c + X_c * beta_anchor;

  Eigen::Map<Eigen::MatrixXd> X_gamma_map(X_gamma.data(), X_gamma.rows(), X_gamma.cols());
  Eigen::Map<Eigen::VectorXd> offset_gamma_map(offset_gamma.data(), offset_gamma.size());

  Rcpp::List gamma_res = opt_pml_hl(
    X_gamma_map,
    offset_gamma_map,
    Y_c,
    fid_c,
    cumsumy_c,
    posind_c,
    posindy_c,
    free_dim,
    X_gamma.rows(),
    fid_c.size() - 1,
    gamma_init,
    sigma_c,
    reml,
    eps,
    ord
  );

  Eigen::VectorXd beta_gamma = Rcpp::as<Eigen::VectorXd>(gamma_res["beta"]);
  Eigen::MatrixXd var_gamma = Rcpp::as<Eigen::MatrixXd>(gamma_res["var"]);

  Eigen::VectorXd beta_full = beta_anchor + K_c * beta_gamma;
  Eigen::MatrixXd var_full = K_c * var_gamma * K_c.transpose();

  if (gamma_res.containsElementNamed("per_subject_stats"))
  {
    Rcpp::List ps = gamma_res["per_subject_stats"];
    if (free_dim > 0)
    {
      Eigen::MatrixXd h_beta_gamma = Rcpp::as<Eigen::MatrixXd>(ps["h_beta"]);
      Eigen::MatrixXd h_eta_beta_gamma = Rcpp::as<Eigen::MatrixXd>(ps["h_eta_beta"]);
      Eigen::MatrixXd h_etaeta_beta_gamma = Rcpp::as<Eigen::MatrixXd>(ps["h_etaeta_beta"]);
      Eigen::MatrixXd h_beta_phi_gamma = Rcpp::as<Eigen::MatrixXd>(ps["h_beta_phi"]);
      Eigen::MatrixXd b_beta_gamma = Rcpp::as<Eigen::MatrixXd>(ps["b_beta"]);
      Eigen::MatrixXd Kmat = K_c;
      ps["h_beta"] = Kmat * h_beta_gamma;
      ps["h_eta_beta"] = Kmat * h_eta_beta_gamma;
      ps["h_etaeta_beta"] = Kmat * h_etaeta_beta_gamma;
      ps["h_beta_phi"] = Kmat * h_beta_phi_gamma;
      ps["b_beta"] = Kmat * b_beta_gamma;

      if (ps.containsElementNamed("h_etaeta_beta_beta"))
      {
        Rcpp::List beta_beta_list = ps["h_etaeta_beta_beta"];
        for (int idx = 0; idx < subject_count; ++idx)
        {
          Eigen::MatrixXd mat_gamma = Rcpp::as<Eigen::MatrixXd>(beta_beta_list[idx]);
          Eigen::MatrixXd mat_full = Kmat * mat_gamma * Kmat.transpose();
          beta_beta_list[idx] = Rcpp::wrap(mat_full);
        }
        ps["h_etaeta_beta_beta"] = beta_beta_list;
      }

      if (ps.containsElementNamed("h_beta_beta"))
      {
        Rcpp::List beta_beta_list = ps["h_beta_beta"];
        for (int idx = 0; idx < subject_count; ++idx)
        {
          Eigen::MatrixXd mat_gamma = Rcpp::as<Eigen::MatrixXd>(beta_beta_list[idx]);
          Eigen::MatrixXd mat_full = Kmat * mat_gamma * Kmat.transpose();
          beta_beta_list[idx] = Rcpp::wrap(mat_full);
        }
        ps["h_beta_beta"] = beta_beta_list;
      }

      if (ps.containsElementNamed("g_beta_beta"))
      {
        Rcpp::List g_beta_beta_list = ps["g_beta_beta"];
        for (int idx = 0; idx < subject_count; ++idx)
        {
          Eigen::MatrixXd mat_gamma = Rcpp::as<Eigen::MatrixXd>(g_beta_beta_list[idx]);
          Eigen::MatrixXd mat_full = Kmat * mat_gamma * Kmat.transpose();
          g_beta_beta_list[idx] = Rcpp::wrap(mat_full);
        }
        ps["g_beta_beta"] = g_beta_beta_list;
      }

      auto transform_vector_list = [&](const char *name)
      {
        if (!ps.containsElementNamed(name))
        {
          return;
        }
        Rcpp::List vec_list = ps[name];
        for (int idx = 0; idx < subject_count; ++idx)
        {
          Eigen::VectorXd vec_gamma = Rcpp::as<Eigen::VectorXd>(vec_list[idx]);
          Eigen::VectorXd vec_full = Kmat * vec_gamma;
          vec_list[idx] = Rcpp::wrap(vec_full);
        }
        ps[name] = vec_list;
      };

      transform_vector_list("h_etaeta_beta_phi");
      transform_vector_list("f_beta_phi");
      transform_vector_list("g_beta_phi");
    }
    else
    {
      const int nb_full = beta_anchor.size();
      ps["h_beta"] = Eigen::MatrixXd::Zero(nb_full, subject_count);
    ps["h_eta_beta"] = Eigen::MatrixXd::Zero(nb_full, subject_count);
    ps["h_etaeta_beta"] = Eigen::MatrixXd::Zero(nb_full, subject_count);
    ps["h_beta_phi"] = Eigen::MatrixXd::Zero(nb_full, subject_count);
    ps["b_beta"] = Eigen::MatrixXd::Zero(nb_full, subject_count);

      auto make_zero_matrix_list = [&](void) {
        Rcpp::List res(subject_count);
        Eigen::MatrixXd zero_mat = Eigen::MatrixXd::Zero(nb_full, nb_full);
        for (int idx = 0; idx < subject_count; ++idx)
        {
          res[idx] = Rcpp::wrap(zero_mat);
        }
        return res;
      };
      ps["h_etaeta_beta_beta"] = make_zero_matrix_list();
      ps["h_beta_beta"] = make_zero_matrix_list();
      ps["g_beta_beta"] = make_zero_matrix_list();

      auto assign_zero_vector_list = [&](const char *name)
      {
        Rcpp::List zero_list(subject_count);
        Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(nb_full);
        for (int idx = 0; idx < subject_count; ++idx)
        {
          zero_list[idx] = Rcpp::wrap(zero_vec);
        }
        ps[name] = zero_list;
      };

      assign_zero_vector_list("h_etaeta_beta_phi");
      assign_zero_vector_list("f_beta_phi");
      assign_zero_vector_list("g_beta_phi");
    }
    if (!ps.containsElementNamed("g_tau_tau"))
    {
      Eigen::VectorXd g_tau_vec = -Rcpp::as<Eigen::VectorXd>(ps["h_etaeta_tau"]);
      ps["g_tau_tau"] = g_tau_vec;
    }
    gamma_res["per_subject_stats"] = ps;
  }

  gamma_res["beta"] = beta_full;
  gamma_res["var"] = var_full;
  return gamma_res;
}

// [[Rcpp::export]]
Rcpp::List opt_pml_nbm(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                   const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                   const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                   const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                   const Eigen::VectorXd & beta_c,
                   const Eigen::VectorXd & sigma_c, const int reml, const double eps, const int ord)
{
  double alpha = sigma_c(0);
  double gamma = sigma_c(1);

  double loglik = 0;

  Eigen::VectorXd logw = Eigen::VectorXd::Constant(k_c,1,0);
  // Eigen::VectorXd logw = logw_c;
  Eigen::VectorXd beta = beta_c;

  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}

  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }

  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  Eigen::VectorXd extb = offset_c + X_c*beta;

  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}

  loglik += logw.dot(cumsumy_c);

  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();

  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}

  loglik = loglik - logw.dot(logw)/alpha/2 - k_c/2*log(alpha);


  // double eps = 1e-6;
  double loglikp = 0;
  double likdif = 0;
  int step = 0;
  Eigen::MatrixXd vb(nb_c,nb_c);
  Eigen::VectorXd vw(k_c);
  Eigen::MatrixXd vb2(nb_c,nb_c);
  Eigen::MatrixXd vwb(k_c,nb_c);
  Eigen::VectorXd gstar_extb_phi(nind_c);
  int stepd = 0;
  int maxstd = 10;
  int maxstep = 50;
  double convd = 0.01;

  //while((step==0)||(likdif>abs(eps*loglik)))
  while((step==0)||((likdif>eps)&&(step<maxstep)))
  {
    step++;

    // double damp = 1;
    Eigen::VectorXd damp_w = Eigen::VectorXd::Constant(k_c,1);
    Eigen::VectorXd damp = Eigen::VectorXd::Constant(nb_c,1);

    gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
    Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
    Eigen::VectorXd dw(k_c);
    for(int i=0;i<k_c;i++)
    {
      dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    dw = cumsumy_c-dw-logw/alpha;

    gstar_extb_phi = gstar_extb_phi.array()/(extb.array()+gamma);

    for(int i=0;i<k_c;i++)
    {
      vw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
    }
    vw = gamma*vw;
    vw = vw.array() + 1/alpha;

    Eigen::MatrixXd xgsetbp = X_c.array().colwise()*gstar_extb_phi.array();

    for(int i=0;i<k_c;i++)
    {
      vwb.row(i) = xgsetbp.block(fid_c[i],0,len(i),nb_c).colwise().sum();
    }
    vwb = gamma*vwb;
    for(int i=0;i<nb_c;i++)
    {
      for(int j=i;j<nb_c;j++)
      {
        vb(i,j) = X_c.col(i).dot(xgsetbp.col(j));
        if(i!=j)
        {
          vb(j,i) = vb(i,j);
        }
      }
    }
    // vb = X_c.transpose()*vb;
    vb = gamma*vb;
    // Eigen::MatrixXd temp = vwb.array().colwise()/vw.array();
    // vb2 = vb - vwb.transpose()*temp;
    Eigen::MatrixXd temp = vwb.array().colwise()/vw.array().sqrt();
    vb2 = vb - temp.transpose()*temp;

    Eigen::VectorXd dwvw = dw.array()/vw.array();
    // Eigen::VectorXd dbvwbvbdw = vb.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd stepbeta = vb2.ldlt().solve(db - vwb.transpose()*dwvw);
    Eigen::VectorXd new_b = beta + stepbeta;
    Eigen::VectorXd steplogw = dwvw - ((vwb*stepbeta).array()/vw.array()).matrix();
    Eigen::VectorXd new_w = logw + steplogw;

    loglikp = loglik;
    loglik = 0;

    extb = offset_c + X_c*new_b;

    for(unsigned int i=0;i<nelem;i++)
    {loglik += extb(posindy_c(i))*Y_c(i);}
    loglik += new_w.dot(cumsumy_c);
    for(int i=0;i<k_c;i++)
    {
      extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
    }
    extb = extb.array().exp();

    extbphil = (extb.array()+gamma).log();
    loglik -= gamma*extbphil.sum();
    for(unsigned int i=0;i<nelem;i++)
    {loglik -= Y_c(i)*extbphil(posindy_c(i));}
    loglik = loglik - new_w.dot(new_w)/alpha/2 - k_c/2*log(alpha);

    likdif = loglik - loglikp;
    stepd = 0;
    double minstep = 40;
    
    while((likdif<0)||(std::isinf(loglik)))
    {
      stepd++;
      minstep = minstep/2;
      
      if(stepd>maxstd)
      {
        likdif = 0;
        loglik = loglikp;
        double mabsdb = db.cwiseAbs().maxCoeff();
        double mabsdw = dw.cwiseAbs().maxCoeff();
        if((mabsdb>convd) || (mabsdw>convd))
        {stepd++;}
        break;
      }
      
      for(int i=0;i<nb_c;i++)
      {
        if((stepbeta(i)<40)&&(stepbeta(i)>-40))
        {
          damp(i)=damp(i)/2;
          new_b(i) = beta(i) + stepbeta(i)*damp(i);
        }else{
          if(stepbeta(i)>0)
          {
            new_b(i) = beta(i) + minstep;
          }else{
            new_b(i) = beta(i) - minstep;
          }
        }
      }
      
      for(int i=0;i<k_c;i++)
      {
        if((steplogw(i)<40)&&(steplogw(i)>-40))
        {
          damp_w(i)=damp_w(i)/2;
          new_w(i) = logw(i) + steplogw(i)*damp_w(i);
        }else{
          if(steplogw(i)>0)
          {
            new_w(i) = logw(i) + minstep;
          }else{
            new_w(i) = logw(i) - minstep;
          }
        }
      }

      loglik = 0;

      extb = offset_c + X_c*new_b;

      for(unsigned int i=0;i<nelem;i++)
      {loglik += extb(posindy_c(i))*Y_c(i);}
      loglik += new_w.dot(cumsumy_c);
      for(int i=0;i<k_c;i++)
      {
        extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+new_w(i);
      }
      extb = extb.array().exp();

      extbphil = (extb.array()+gamma).log();
      loglik -= gamma*extbphil.sum();
      for(unsigned int i=0;i<nelem;i++)
      {loglik -= Y_c(i)*extbphil(posindy_c(i));}
      loglik = loglik - new_w.dot(new_w)/alpha/2 - k_c/2*log(alpha);

      likdif = loglik - loglikp;
    }
    
    beta = new_b;
    logw = new_w;
  }

  double logdet;
  if(reml==1)
  {
    logdet = vw.array().abs().log().sum();
    double absdet = vb2.determinant();
    if(absdet<0)
      absdet *= -1;
    logdet += log(absdet);
  }else{
    logdet = vw.array().abs().log().sum();
  }
  // double logdet = vw.array().log().sum();
  
  double sec_ord = 0;

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("logw") = logw,
                            //Rcpp::Named("var") = vb2.inverse(),
                            Rcpp::Named("var") = vb2,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("loglikp") = loglikp,
                            Rcpp::Named("logdet") = logdet,
                            Rcpp::Named("iter") = step,
                            Rcpp::Named("damp") = stepd,
                            Rcpp::Named("second") = sec_ord);

}


// [[Rcpp::export]]
Rcpp::List pml_ll_der_eigen(const Eigen::Map<Eigen::MatrixXd> & X_c, const Eigen::Map<Eigen::VectorXd> & offset_c,
                            const Eigen::VectorXd & Y_c, const Eigen::VectorXi & fid_c,
                            const Eigen::VectorXd & cumsumy_c,const Eigen::VectorXi & posind_c,
                            const Eigen::VectorXi & posindy_c, const int nb_c, const int nind_c, const int k_c,
                            const Eigen::VectorXd & beta_c, const Eigen::VectorXd & logw_c, const Eigen::VectorXd & sigma_c)
{
  double alpha = sigma_c(0);
  double gamma = sigma_c(1);

  double loglik = 0;

  Eigen::VectorXd logw = logw_c;
  Eigen::VectorXd beta = beta_c;

  Eigen::VectorXd gstar = Eigen::VectorXd::Constant(nind_c,1,gamma);
  unsigned int nelem = posindy_c.size();
  for(unsigned int i=0;i<nelem;i++)
  {gstar(posindy_c(i)) += Y_c(i);}

  Eigen::ArrayXd len(k_c);
  for(int i=0;i<k_c;i++)
  {
    len(i) = fid_c[i+1]-fid_c[i];
  }

  Eigen::VectorXd yx = Eigen::VectorXd::Zero(nb_c);
  for(unsigned int i=0;i<nelem;i++)
  {yx += X_c.row(posindy_c(i)).transpose()*Y_c(i);}

  Eigen::VectorXd extb = offset_c + X_c*beta;

  Eigen::VectorXd w = logw.array().exp();

  for(unsigned int i=0;i<nelem;i++)
  {loglik += extb(posindy_c(i))*Y_c(i);}

  loglik += logw.dot(cumsumy_c);

  for(int i=0;i<k_c;i++)
  {
    extb.segment(fid_c(i),len(i)) = extb.segment(fid_c(i),len(i)).array()+logw(i);
  }
  extb = extb.array().exp();

  Eigen::ArrayXd extbphil = (extb.array()+gamma).log();
  loglik -= gamma*extbphil.sum();
  for(unsigned int i=0;i<nelem;i++)
  {loglik -= Y_c(i)*extbphil(posindy_c(i));}

  loglik += alpha*logw.sum() - alpha*w.sum();

  Eigen::VectorXd gstar_extb_phi = gstar.array()/(1+gamma/extb.array());
  Eigen::VectorXd db = yx - X_c.transpose()*gstar_extb_phi;
  Eigen::VectorXd dw(k_c);
  for(int i=0;i<k_c;i++)
  {
    dw(i) = gstar_extb_phi.segment(fid_c(i),len(i)).sum();
  }
  dw = cumsumy_c - dw - alpha*w;
  dw = dw.array() + alpha;

  Eigen::VectorXd gr(k_c+nb_c);
  gr.head(nb_c) = db;
  gr.tail(k_c) = dw;

  return Rcpp::List::create(Rcpp::Named("fn") = -loglik,
                            Rcpp::Named("gr") = -gr);

}


// [[Rcpp::export]]
Rcpp::List ptmg_ll_der_hes_eigen_per_subject(const Eigen::Map<Eigen::MatrixXd> &X_c, const Eigen::Map<Eigen::VectorXd> &offset_c,
                                             const Eigen::VectorXd &Y_c, const Eigen::VectorXi &fid_c,
                                             const Eigen::VectorXd &cumsumy_c, const Eigen::VectorXi &posind_c,
                                             const Eigen::VectorXi &posindy_c, const int nb_c, const int nind_c, const int k_c,
                                             const Eigen::VectorXd &beta_c, const Eigen::VectorXd &sigma_c)
{
  double exps = exp(sigma_c(0));
  double exps_m = exps - 1;
  double exps_s = sqrt(exps);
  double alpha = 1 / exps_m;
  double lambda = alpha / exps_s;
  double gamma = sigma_c(1);
  double log_lambda = log(lambda);
  double log_gamma = log(gamma);

  exps_m = pow(exps_m, 2);
  double alpha_pr = -exps / exps_m;
  double lambda_pr = (1 - 3 * exps) / (2 * exps_s * exps_m);
  double alpha_dpr = 2 * exps * exps / (exps_m * (exps - 1)) - exps / exps_m;
  double lambda_dpr = -3 * exps / (2 * exps_s * exps_m) + (3 * exps - 1) / (4 * exps_s * exps_m) + (3 * exps - 1) * exps_s / (exps_m * (exps - 1));

  double term1 = 0;
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nb_c + 2, nb_c + 2);
  double hes_sigma = alpha_dpr * log_lambda + 2 * alpha_pr * lambda_pr / lambda - alpha * pow(lambda_pr / lambda, 2) + alpha / lambda * lambda_dpr;
  hes_sigma = hes_sigma * k_c;

  Eigen::VectorXd extb = offset_c + X_c * beta_c;
  unsigned int nelem = posindy_c.size();
  for (unsigned int i = 0; i < nelem; i++)
  {
    term1 += extb(posindy_c(i)) * Y_c(i);
  }

  extb = extb.array().exp();
  Eigen::ArrayXd cumsumxtb(k_c);
  Eigen::ArrayXd len(k_c);
  for (int i = 0; i < k_c; i++)
  {
    len(i) = fid_c[i + 1] - fid_c[i];
    cumsumxtb(i) = extb.segment(fid_c(i), len(i)).sum();
  }

  Eigen::VectorXd ystar = cumsumy_c.array() + alpha;

  Eigen::ArrayXd mustar = cumsumxtb + lambda;
  Eigen::VectorXd mustar_log = mustar.log().matrix();

  term1 -= ystar.dot(mustar_log);

  term1 += k_c * alpha * log_lambda;
  term1 += nind_c * gamma * log_gamma;

  Eigen::VectorXd sum_elgcp(nind_c);
  Eigen::ArrayXd ymustar = ystar.array() / mustar;
  Eigen::ArrayXd ymumustar = ymustar / mustar;
  Eigen::VectorXd gstar_phiymustar = Eigen::VectorXd::Constant(nind_c, 1, gamma);
  for (unsigned int i = 0; i < nelem; i++)
  {
    gstar_phiymustar(posindy_c(i)) += Y_c(i);
  }

  for (int i = 0; i < k_c; i++)
  {
    sum_elgcp.segment(fid_c(i), len(i)) = ymustar(i) * extb.segment(fid_c(i), len(i));
  }

  term1 += sum_elgcp.sum();

  sum_elgcp = sum_elgcp.array() + gamma;
  Eigen::ArrayXd tempa = 1 / sum_elgcp.array();
  hessian(nb_c + 1, nb_c + 1) = nind_c / gamma - 2 * tempa.sum();
  gstar_phiymustar = gstar_phiymustar.array() * tempa;

  sum_elgcp = sum_elgcp.array().log();
  double slpey = sum_elgcp.sum();
  term1 -= gamma * slpey;
  for (unsigned int i = 0; i < nelem; i++)
  {
    term1 -= Y_c(i) * sum_elgcp(posindy_c(i));
  }
  term1 = term1 * (-1);

  Eigen::VectorXd der = Eigen::VectorXd::Zero(nb_c + 2);

  Eigen::MatrixXd xexb = X_c.array().colwise() * extb.array();
  Eigen::MatrixXd xexb_f(k_c, nb_c);
  Eigen::MatrixXd dbeta_41(k_c, nb_c);

  // Initialize obs_to_subject as Eigen::VectorXi
  Eigen::VectorXi obs_to_subject(X_c.rows());
  obs_to_subject.setConstant(-1); // Initialize all elements to -1

  for (int subject = 0; subject < k_c; subject++)
  {
    int start = fid_c[subject];
    int end = fid_c[subject + 1];
    int lent = end - start;

    // Aggregate sums over each subject/group
    xexb_f.row(subject) = xexb.block(start, 0, lent, nb_c).colwise().sum();
    dbeta_41.row(subject) = gstar_phiymustar.segment(start, lent).transpose() * xexb.block(start, 0, lent, nb_c);

    // Assign subject index to obs_to_subject
    for (int pos = start; pos < end; pos++)
    { // pos < end to include 'start' and exclude 'end'
      obs_to_subject(pos) = subject;
    }
  }
  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();

  Eigen::ArrayXd dbeta_42(k_c);
  for (int i = 0; i < k_c; i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    // dbeta_41.row(i) = gstar_phiymustar.segment(start,lent).transpose()*xexb.block(start,0,lent,nb_c);
    dbeta_42[i] = gstar_phiymustar.segment(start, lent).dot(extb.segment(start, lent));
    // dbeta_42[i] = gstar_phiymustar.segment(start,lent).sum();
  }

  dbeta_42 = dbeta_42 - cumsumxtb;

  Eigen::ArrayXd ymumustar_dbeta_csxtb = ymumustar * dbeta_42;
  Eigen::VectorXd db = xexb_f * ymumustar_dbeta_csxtb.matrix() - dbeta_41 * ymustar.matrix();
  for (unsigned int i = 0; i < nelem; i++)
  {
    db += X_c.row(posindy_c(i)).transpose() * Y_c(i);
  }

  // Initialize the per-subject gradient matrix
  Eigen::MatrixXd db_sub = Eigen::MatrixXd::Zero(nb_c, k_c);

  // Ensure ymumustar_dbeta_csxtb and ymustar are column vectors
  Eigen::VectorXd ymumustar_dbeta_csxtb_vec = ymumustar_dbeta_csxtb.matrix(); // (k_c, 1)
  Eigen::VectorXd ymustar_vec = ymustar.matrix();                             // (k_c, 1)

  // Perform element-wise multiplications per column
  for (int k = 0; k < k_c; ++k)
  {
    db_sub.col(k) = xexb_f.col(k) * ymumustar_dbeta_csxtb_vec(k) - dbeta_41.col(k) * ymustar_vec(k);
  }

  // Accumulate contributions from positive observations
  for (unsigned int i = 0; i < nelem; i++)
  {
    int pos = posindy_c(i);            // Observation index in X_c
    int subject = obs_to_subject(pos); // Retrieve the subject index
    // Update the per-subject gradient
    db_sub.col(subject) += X_c.row(pos).transpose() * Y_c(i);
  }

  // // Compute the column-wise sum of db_sub
  // Eigen::VectorXd db_sum = db_sub.rowwise().sum();

  // // Print specific elements for manual inspection
  // std::cout << "db_sum (Sum of db_sub):\n" << db_sum << std::endl;
  // std::cout << "db (Aggregate Gradient):\n" << db << std::endl;

  double dtau = 0;
  double ldm = log_lambda * k_c - mustar_log.sum();
  double adlmy = exps_s * k_c - ymustar.sum();

  Eigen::ArrayXd imustar = 1 / mustar;
  Eigen::VectorXd hes_b_l = dbeta_42 * imustar;
  double dbim = hes_b_l.sum();
  dtau = dtau - alpha_pr * dbim;
  dtau = dtau + lambda_pr * ymumustar_dbeta_csxtb.sum();
  dtau += alpha_pr * ldm + lambda_pr * adlmy;

  // not necessary
  hes_sigma -= dbim * alpha_dpr - 2 * alpha_pr * lambda_pr * hes_b_l.dot(imustar.matrix()) - lambda_dpr * ymumustar_dbeta_csxtb.sum() + 2 * lambda_pr * lambda_pr * (ymumustar_dbeta_csxtb * imustar).sum();

  double dtau2 = log_gamma * nind_c + nind_c - slpey - gstar_phiymustar.sum();

  // // Initialize per-subject dtau vector
  // Eigen::VectorXd dtau_sub = Eigen::VectorXd::Zero(k_c); // Size (k_c)

  // Compute ldm_sub and adlmy_sub as Eigen::VectorXd
  Eigen::VectorXd ldm_sub = Eigen::VectorXd::Constant(k_c, log_lambda).array() - mustar_log.array();
  Eigen::VectorXd adlmy_sub = Eigen::VectorXd::Constant(k_c, exps_s).array() - ymustar.array();

  // // Compute imustar (inverse of mustar)
  // Eigen::ArrayXd imustar = 1/mustar;
  // Eigen::VectorXd hes_b_l = dbeta_42*imustar;      // hes_b_l[i] = dbeta_42[i] * imustar[i]

  // Compute dtau_sub as Eigen::VectorXd
  Eigen::VectorXd dtau_sub = (-alpha_pr * hes_b_l) + (lambda_pr * ymumustar_dbeta_csxtb.matrix()) + (alpha_pr * ldm_sub) + (lambda_pr * adlmy_sub);

  // // Compute the column-wise sum of dtau_sub
  // double dtau_sum = dtau_sub.sum();

  // // Print specific elements for manual inspection
  // std::cout << "dtau_sum (Sum of dtau_sum):\n" << dtau_sum << std::endl;
  // std::cout << "dtau (Aggregate Gradient):\n" << dtau << std::endl;

  // Compute per-subject dtau2_sub as Eigen::VectorXd
  Eigen::VectorXd dtau2_sub = Eigen::VectorXd::Zero(k_c);

  for (int i = 0; i < k_c; ++i)
  {
    int start = fid_c[i];
    int lent = len[i];

    double n_i = static_cast<double>(lent);

    double slpey_sub_i = sum_elgcp.segment(start, lent).sum();
    double gstar_phiymustar_sub_i = gstar_phiymustar.segment(start, lent).sum();

    dtau2_sub[i] = (log_gamma * n_i) + n_i - slpey_sub_i - gstar_phiymustar_sub_i;
  }

  // Compute aggregate dtau2
  // double dtau2_sum = dtau2_sub.sum();

  // // Print specific elements for manual inspection
  // std::cout << "dtau2_sum (Sum of dtau2_sub):\n" << dtau2_sum << std::endl;
  // std::cout << "dtau2 (Aggregate Gradient):\n" << dtau2 << std::endl;

  // Eigen::VectorXd hes_sigma_beta = -alpha_pr*(xexb_f*imustar.matrix()) + lambda_pr*(xexb_f*ymumustar.matrix());
  Eigen::VectorXd hes_sigma_beta = -alpha_pr * (dbeta_41 * imustar.matrix()) + lambda_pr * (dbeta_41 * ymumustar.matrix());
  hes_sigma_beta -= -alpha_pr * (xexb_f * hes_b_l.cwiseProduct(imustar.matrix())) + 2 * lambda_pr * (xexb_f * hes_b_l.cwiseProduct(ymumustar.matrix()));

  hes_sigma = hes_sigma - (2 * alpha_pr * lambda_pr * imustar.sum() + alpha_dpr * mustar_log.sum() + lambda_dpr * ymustar.sum() - lambda_pr * lambda_pr * ymumustar.sum());
  Eigen::VectorXd apmmlpymm = alpha_pr * imustar - lambda_pr * ymumustar;

  hes_b_l = gstar_phiymustar.array() * tempa;
  hessian(nb_c + 1, nb_c + 1) = hessian(nb_c + 1, nb_c + 1) + hes_b_l.sum();

  Eigen::VectorXd hes_b_l_f(k_c);
  Eigen::MatrixXd hes_b_l_xexb = xexb.array().colwise() * (hes_b_l.array() - tempa);
  Eigen::MatrixXd hes_b_l_xexb_f(k_c, nb_c);
  for (int i = 0; i < k_c; i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start, 0, lent, nb_c).colwise().sum();
    // hes_b_l_f(i) = hes_b_l.segment(start,lent).sum();
  }
  Eigen::VectorXd hes_gamma_beta = hes_b_l_xexb_f.transpose() * ymustar.matrix();

  hes_b_l = hes_b_l.array() * extb.array();
  Eigen::VectorXd hes_b2_l = hes_b_l - extb.cwiseProduct(tempa.matrix());
  for (int i = 0; i < k_c; i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_f(i) = hes_b2_l.segment(start, lent).sum();
  }
  hes_gamma_beta = hes_gamma_beta - xexb_f * hes_b_l_f.cwiseProduct(ymumustar.matrix());
  hessian.row(nb_c + 1).head(nb_c) = hes_gamma_beta;
  hessian.col(nb_c + 1).head(nb_c) = hes_gamma_beta;
  hessian(nb_c + 1, nb_c) = apmmlpymm.dot(hes_b_l_f);
  hessian(nb_c, nb_c + 1) = hessian(nb_c + 1, nb_c);

  hes_b2_l = hes_b_l.array() * extb.array();

  Eigen::VectorXd tempb;
  Eigen::VectorXd tempc;
  hes_b_l_xexb = xexb.array().colwise() * hes_b_l.array();
  Eigen::VectorXd hes_b2_l_f(k_c);
  for (int i = 0; i < k_c; i++)
  {
    int start = fid_c[i];
    int lent = len(i);
    hes_b_l_xexb_f.row(i) = hes_b_l_xexb.block(start, 0, lent, nb_c).colwise().sum();
    hes_b2_l_f(i) = hes_b2_l.segment(start, lent).sum();
  }
  // hes_b_l_xexb_f.transposeInPlace();

  hes_sigma_beta += hes_b_l_xexb_f.transpose() * (apmmlpymm.cwiseProduct(ymustar.matrix())) - xexb_f * (apmmlpymm.cwiseProduct(hes_b2_l_f.cwiseProduct(ymumustar.matrix())));
  hessian.row(nb_c).head(nb_c) = hes_sigma_beta;
  hessian.col(nb_c).head(nb_c) = hes_sigma_beta;

  apmmlpymm = apmmlpymm.cwiseProduct(apmmlpymm);
  hes_sigma += hes_b2_l_f.dot(apmmlpymm);
  hessian(nb_c, nb_c) = hes_sigma;

  xexb_f.transposeInPlace();
  dbeta_41.transposeInPlace();
  apmmlpymm = ymustar.cwiseProduct(ymustar);
  Eigen::VectorXd ymuymuimu = apmmlpymm.cwiseProduct(imustar.matrix());
  hes_b2_l_f = hes_b2_l_f.array() * ymuymuimu.array() * imustar;

  for (int i = 0; i < nb_c; i++)
  {
    for (int j = i; j < nb_c; j++)
    {
      tempb = xexb.col(i).cwiseProduct(X_c.col(j));
      tempc = tempb.cwiseProduct(hes_b_l);

      for (int k = 0; k < k_c; k++)
      {
        int start = fid_c[k];
        int lent = len(k);

        hessian(i, j) -= tempb.segment(start, lent).dot(gstar_phiymustar.segment(start, lent)) * ymustar(k);
        hessian(i, j) += tempc.segment(start, lent).sum() * apmmlpymm(k);
        // not necessary
        hessian(i, j) += ymumustar_dbeta_csxtb(k) * tempb.segment(start, lent).sum();
      }
      hes_b_l_f = dbeta_41.col(i).cwiseProduct(xexb_f.col(j)) + dbeta_41.col(j).cwiseProduct(xexb_f.col(i));
      hessian(i, j) += hes_b_l_f.dot(ymumustar.matrix());
      hessian(i, j) -= ymuymuimu.dot(hes_b_l_xexb_f.col(i).cwiseProduct(xexb_f.col(j)) + hes_b_l_xexb_f.col(j).cwiseProduct(xexb_f.col(i)));
      tempb = xexb_f.col(i).cwiseProduct(xexb_f.col(j));
      hessian(i, j) -= tempb.dot(ymumustar.matrix());
      hessian(i, j) += hes_b2_l_f.dot(tempb);
      // not necessary
      hessian(i, j) += (-2) * tempb.dot((ymumustar_dbeta_csxtb * imustar).matrix());

      if (i != j)
      {
        hessian(j, i) = hessian(i, j);
      }
    }
  }

  der.segment(0, nb_c) = db;
  der[nb_c] = dtau;
  der[nb_c + 1] = dtau2;

  // Combine db_sub, dtau_sub, and dtau2_sub into per-subject gradients
  Eigen::MatrixXd per_subject_gradients(nb_c + 2, k_c); // Size: (nb_c + 2) x k_c

  for (int k = 0; k < k_c; ++k)
  {
    // The first nb_c elements are db_sub.col(k)
    per_subject_gradients.col(k).head(nb_c) = db_sub.col(k);
    // The next element is dtau_sub(k)
    per_subject_gradients(nb_c, k) = dtau_sub(k);
    // The last element is dtau2_sub(k)
    per_subject_gradients(nb_c + 1, k) = dtau2_sub(k);
  }

  // Compute cumsumy_alpha
  Eigen::VectorXd cumsumy_alpha = cumsumy_c.array() + alpha;

  // Compute Y_plus_gamma
  Eigen::VectorXd Y_plus_gamma = Y_c.array() + gamma;

  // Modify the return statement to include per_subject_gradients
  return Rcpp::List::create(
      Rcpp::Named("fn") = term1,
      Rcpp::Named("gr") = (-1) * der,
      Rcpp::Named("hes") = (-1) * hessian,
      Rcpp::Named("per_subject_gradients") = (-1) * per_subject_gradients,
      Rcpp::Named("cumsumy_alpha") = cumsumy_alpha,
      Rcpp::Named("Y_plus_gamma") = Y_plus_gamma);
}
