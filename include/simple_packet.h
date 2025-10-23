/**
 * @file simple_packet.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a struct for unformed SCHP5 packets
 * @version 0.1
 * @date 2025-10-22
 */
# ifndef SIMPLE_PACKET_H
# define SIMPLE_PACKET_H

#include <cstdint>
#include <vector>

/**
 * @brief Implements a packet compliant to the SCHP5 protocol
 */
struct SimpPacket {

    /**
     * @brief Describes the packet type
     */
    uint8_t controlCode;

    /**
     * @brief Describes how many bytes are contained in the packet
     */
    uint8_t dataLength;

    /**
     * @brief The address of the sender of this packet
     */
    uint8_t recieverAddress;

    /**
     * @brief The address of the sender of this packet
     */
    uint8_t senderAddress;


    /**
     * @brief The sequence number of the current packet
     */
    uint8_t sequenceNumber;

    /** 
     * @brief The data content of this packet
    */
   std::vector<uint8_t> data;

   /**
    * @brief True checksum returned good 
    */
   bool erc;

};

#endif