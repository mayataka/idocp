#include "idocp/hybrid/ocp_discretizer.hpp"

#include <iostream>

namespace idocp {

void OCPDiscretizer::showInfo() const {
  std::cout << "----- The discretized optimal control problem (OCP) -----" << std::endl;
  std::cout << "T = " << T_ << std::endl;
  std::cout << "N_ideal = " << N_ideal() << std::endl;
  std::cout << "N = " << N() << std::endl;
  std::cout << "N_impulse = " << N_impulse() << std::endl;
  std::cout << "N_lift = " << N_lift() << std::endl;
  std::cout << "N_all = " << N_all() << std::endl;
  for (int i=0; i<N(); ++i) {
    std::cout << "i = " << i << ": t = " << t(i) << ": dt = " << dt(i) << std::endl; 
    if (isTimeStageBeforeImpulse(i)) {
      const int index = impulseIndexAfterTimeStage(i);
      std::cout << "impulse = " << index 
                << ": t =  " << t_impulse(index) 
                << ": dt_aux = " << dt_aux(index) << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int index = liftIndexAfterTimeStage(i);
      std::cout << "lift = " << index 
                << ": t =  " << t_lift(index) 
                << ": dt_lift = " << dt_lift(index) << std::endl;
    }
  }
  std::cout << "i = " << N() << ": t = " << t(N()) << std::endl; 
}


void OCPDiscretizer::showInfo(const ContactSequence& contact_sequence) const {
  std::cout << "----- The discretized optimal control problem (OCP) -----" << std::endl;
  std::cout << "T = " << T_ << std::endl;
  std::cout << "N_ideal = " << N_ideal() << std::endl;
  std::cout << "N = " << N() << std::endl;
  std::cout << "N_impulse = " << N_impulse() << std::endl;
  std::cout << "N_lift = " << N_lift() << std::endl;
  std::cout << "N_all = " << N_all() << std::endl;
  for (int i=0; i<N(); ++i) {
    std::cout << "i = " << i << ": t = " << t(i) << ": dt = " << dt(i) 
              << ": active contacts = [";
    const auto contact_status = contact_sequence.contactStatus(contactPhase(i)); 
    for (int j=0; j<contact_status.maxPointContacts(); ++j) {
      if (contact_status.isContactActive(j)) {
        std::cout << j << " ";
      }
    }
    std::cout << "]" << std::endl;
    // std::cout << "contact points = [";
    // for (int j=0; j<contact_status.maxPointContacts(); ++j) {
    //   if (contact_status.isContactActive(j)) {
    //     std::cout << contact_status.contactPoint(j).transpose() << "; ";
    //   }
    // }
    // std::cout << "]" << std::endl;
    if (isTimeStageBeforeImpulse(i)) {
      const int index = impulseIndexAfterTimeStage(i);
      std::cout << "impulse = " << index 
                << ": t =  " << t_impulse(index) 
                << ": dt_aux = " << dt_aux(index)
                << ": active impulse = [";
      const auto impulse_status = contact_sequence.impulseStatus(index); 
      for (int j=0; j<impulse_status.maxPointContacts(); ++j) {
        if (impulse_status.isImpulseActive(j)) {
          std::cout << j << " ";
        }
      }
      std::cout << "]" << std::endl;
      // std::cout << "impulse points = [";
      // for (int j=0; j<impulse_status.maxPointContacts(); ++j) {
      //   if (impulse_status.isImpulseActive(j)) {
      //     std::cout << impulse_status.contactPoint(j).transpose() << "; ";
      //   }
      // }
      // std::cout << "]" << std::endl;
      std::cout << ": active contacts aux = [";
      const auto contact_status = contact_sequence.contactStatus(contactPhaseAfterImpulse(index)); 
      for (int j=0; j<contact_status.maxPointContacts(); ++j) {
        if (contact_status.isContactActive(j)) {
          std::cout << j << " ";
        }
      }
      std::cout << "]" << std::endl;
      std::cout << "contact points = [";
      for (int j=0; j<contact_status.maxPointContacts(); ++j) {
      if (contact_status.isContactActive(j)) {
        std::cout << contact_status.contactPoint(j).transpose() << "; ";
      }
    }
    }
    else if (isTimeStageBeforeLift(i)) {
      const int index = liftIndexAfterTimeStage(i);
      std::cout << "lift = " << index 
                << ": t =  " << t_lift(index) 
                << ": dt_lift = " << dt_lift(index)
                << ": active contacts lift = [";
      const auto contact_status = contact_sequence.contactStatus(contactPhaseAfterLift(index)); 
      for (int j=0; j<contact_status.maxPointContacts(); ++j) {
        if (contact_status.isContactActive(j)) {
          std::cout << j << " ";
        }
      }
      std::cout << "]" << std::endl;
    }
  }
  std::cout << "i = " << N() << ": t = " << t(N()) << std::endl; 
}

} // namespace idocp 