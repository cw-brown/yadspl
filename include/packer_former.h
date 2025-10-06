/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for forming SCHP5 packets
 * @version 0.1
 * @date 2025-06-18
 */
# ifndef PACKET_FORMER_H
# define PACKET_FORMER_H

#include "ring.hpp"

/**
 * @brief Forms packets according to the SCHP5 protocol for packet radios
 */
class PacketFormer {

    private:

    /**
     * @brief Pointer to the buffer where the packet former reads in arbitrary binary data formatted as a buffer of bytes
     */
    std::ring<char> * inputBuffer;

    public:

    /**
     * @brief Construct a new packet former
     */
    PacketFormer(std::ring<char> * inputBuffer);

};

/**
 * @brief Implements a packet compliant to the SCHP5 protocol
 */
class Packet {

    private:

    // Fields
    /**
     * @brief The flag at the start and stop of each packet
     */
    static const char flag = 0b01111110;

    /**
     * @brief Describes the packet type
     */
    char controlCode;

    /**
     * @brief Describes how many bytes are contained in the packet
     */
    char dataLength;

    /**
     * @brief The adress of the sender of this packet
     */
    

};

#endif