/**
 * @file preceiver.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for receiving SCHP5 Packets
 * @version 0.1
 * @date 2025-11-04
 */

 # ifndef PRECEIVER_H
 # define PRECEIVER_H

#include <cstdint>
#include "packet.h"
#include "simple_packet.h"
#include "packet_former.h"
#include "packet_reformer.h"
#include "data_extractor.h"

// State definitions
enum State {
    START,
    WAIT_DATA,
    UNPACKET,
    SEND,
    WAIT_ACK,
};

/** 
 * @brief Handles control flow and packeting for SCHP5 data sink
 */
class PacketReceiver {

};
 
 # endif