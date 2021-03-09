#ifndef IDOCP_IMPULSE_STATUS_HXX_
#define IDOCP_IMPULSE_STATUS_HXX_

#include "idocp/robot/impulse_status.hpp"

#include <cassert>
#include <random>


namespace idocp {

inline ImpulseStatus::ImpulseStatus(const int max_point_contacts)
  : impulse_status_(max_point_contacts) {
}


inline ImpulseStatus::ImpulseStatus() 
  : impulse_status_() {
}
 

inline ImpulseStatus::~ImpulseStatus() {
}


inline bool ImpulseStatus::operator==(const ImpulseStatus& other) const {
  assert(other.maxPointContacts() == maxPointContacts());
  for (int i=0; i<maxPointContacts(); ++i) {
    if (other.isImpulseActive(i) != isImpulseActive(i)) {
      return false;
    }
    if (!other.contactPoints()[i].isApprox(contactPoints()[i])) {
      return false;
    }
  }
  return true;
}


inline bool ImpulseStatus::operator!=(const ImpulseStatus& other) const {
  return !(*this == other);
}


inline bool ImpulseStatus::isImpulseActive(const int contact_index) const {
  return impulse_status_.isContactActive(contact_index);
}


inline const std::vector<bool>& ImpulseStatus::isImpulseActive() const {
  return impulse_status_.isContactActive();
}


inline bool ImpulseStatus::hasActiveImpulse() const {
  return impulse_status_.hasActiveContacts();
}


inline int ImpulseStatus::dimf() const {
  return impulse_status_.dimf();
}


inline int ImpulseStatus::maxPointContacts() const {
  return impulse_status_.maxPointContacts();
}


inline void ImpulseStatus::setActivity(
    const ContactStatus& pre_contact_status, 
    const ContactStatus& post_contact_status) {
  assert(pre_contact_status.maxPointContacts() == maxPointContacts());
  assert(post_contact_status.maxPointContacts() == maxPointContacts());
  for (int i=0; i<impulse_status_.maxPointContacts(); ++i) {
    if (pre_contact_status.isContactActive(i)) {
      deactivateImpulse(i);
    }
    else {
      if (post_contact_status.isContactActive(i)) {
        activateImpulse(i);
      }
      else {
        deactivateImpulse(i);
      }
    }
  }
}


inline void ImpulseStatus::setActivity(
    const std::vector<bool>& is_impulse_active) {
  impulse_status_.setActivity(is_impulse_active);
}


inline void ImpulseStatus::activateImpulse(const int contact_index) {
  impulse_status_.activateContact(contact_index);
}


inline void ImpulseStatus::deactivateImpulse(const int contact_index) {
  impulse_status_.deactivateContact(contact_index);
}


inline void ImpulseStatus::activateImpulse(
    const std::vector<int>& contact_indices) {
  impulse_status_.activateContacts(contact_indices);
} 


inline void ImpulseStatus::activateImpulse() {
  impulse_status_.activateContacts();
} 


inline void ImpulseStatus::deactivateImpulse(
    const std::vector<int>& contact_indices) {
  impulse_status_.deactivateContacts(contact_indices);
}


inline void ImpulseStatus::deactivateImpulse() {
  impulse_status_.deactivateContacts();
}


inline void ImpulseStatus::setContactPoint(
    const int contact_index, const Eigen::Vector3d& contact_point) {
  impulse_status_.setContactPoint(contact_index, contact_point);
}


inline void ImpulseStatus::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  impulse_status_.setContactPoints(contact_points);
}


inline const std::vector<Eigen::Vector3d>& 
ImpulseStatus::contactPoints() const {
  return impulse_status_.contactPoints();
}


inline const Eigen::Vector3d& ImpulseStatus::contactPoint(const int j) const {
  return impulse_status_.contactPoint(j);
}


inline void ImpulseStatus::setRandom() {
  impulse_status_.setRandom();
}

} // namespace idocp

#endif // IDOCP_IMPULSE_STATUS_HXX_ 